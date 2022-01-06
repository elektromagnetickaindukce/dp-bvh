#ifndef BVH_UTILITIES_HPP
#define BVH_UTILITIES_HPP

#include <cstring>
#include <cstdint>
#include <atomic>
#include <memory>
#include <queue>
#include <algorithm>
#include <cmath>
#include <climits>
#include <math.h>
#include <type_traits>
#include <iostream>

#include "bvh/bounding_box.hpp"

namespace bvh {

	/// Safe function to reinterpret the bits of the given value as another type.
	template <typename To, typename From>
	To as(From from) {
		static_assert(sizeof(To) == sizeof(From));
		To to;
		std::memcpy(&to, &from, sizeof(from));
		return to;
	}

	/// Equivalent to copysign(x, x * y).
	inline float product_sign(float x, float y) {
		return as<float>(as<uint32_t>(x) ^ (as<uint32_t>(y) & UINT32_C(0x80000000)));
	}

	/// Equivalent to copysign(x, x * y).
	inline double product_sign(double x, double y) {
		return as<double>(as<uint64_t>(x) ^ (as<uint64_t>(y) & UINT64_C(0x8000000000000000)));
	}

	inline float fast_multiply_add(float x, float y, float z) {
#ifdef FP_FAST_FMAF
		return std::fmaf(x, y, z);
#else
		return x * y + z;
#endif
	}

	inline double fast_multiply_add(double x, double y, double z) {
#ifdef FP_FAST_FMA
		return std::fma(x, y, z);
#else
		return x * y + z;
#endif
	}

	/// Returns the mininum of two values.
	/// Guaranteed to return a non-NaN value if the right hand side is not a NaN.
	template <typename T>
	const T& robust_min(const T& x, const T& y) {
		return x < y ? x : y;
	}

	/// Returns the maximum of two values.
	/// Guaranteed to return a non-NaN value if the right hand side is not a NaN.
	template <typename T>
	const T& robust_max(const T& x, const T& y) {
		return x > y ? x : y;
	}

	template <typename T>
	void atomic_max(std::atomic<T>& x, T y) {
		auto z = x.load();
		while (z < y && !x.compare_exchange_weak(z, y));
	}

	/// Templates that contains signed and unsigned integer types of the given number of bits.
	template <size_t Bits>
	struct SizedIntegerType {
		static_assert(Bits <= 8);
		using Signed = int8_t;
		using Unsigned = uint8_t;
	};

	template <>
	struct SizedIntegerType<64> {
		using Signed = int64_t;
		using Unsigned = uint64_t;
	};

	template <>
	struct SizedIntegerType<32> {
		using Signed = int32_t;
		using Unsigned = uint32_t;
	};

	template <>
	struct SizedIntegerType<16> {
		using Signed = int16_t;
		using Unsigned = uint16_t;
	};

	/// Adds the given number of ULPs (Unit in the Last Place) to the floating-point argument.
	template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
	T add_ulp_magnitude(T x, unsigned ulps) {
		using U = typename SizedIntegerType<sizeof(T) * CHAR_BIT>::Unsigned;
		return std::isfinite(x) ? as<T>(as<U>(x) + ulps) : x;
	}

	/// Computes the (rounded-up) compile-time log in base-2 of an unsigned integer.
	inline constexpr size_t round_up_log2(size_t i, size_t p = 0) {
		return (size_t(1) << p) >= i ? p : round_up_log2(i, p + 1);
	}

	/// Returns the number of bits that are equal to zero,
	/// starting from the most significant one.
	template <typename T, std::enable_if_t<std::is_unsigned<T>::value, int> = 0>
	size_t count_leading_zeros(T value) {
		static constexpr size_t bit_count = sizeof(T) * CHAR_BIT;
		size_t a = 0;
		size_t b = bit_count;
		auto all = T(-1);
		for (size_t i = 0; i < round_up_log2(bit_count); i++) {
			auto m = (a + b) / 2;
			auto mask = all << m;
			if (value & mask) a = m + 1;
			else              b = m;
		}
		return bit_count - b;
	}

	/// Shuffles primitives such that the primitive at index i is `primitives[indices[i]]`.
	template <typename Primitive>
	std::unique_ptr<Primitive[]> shuffle_primitives(const Primitive* primitives, const size_t* indices, size_t primitive_count) {
		auto primitives_copy = std::make_unique<Primitive[]>(primitive_count);
#pragma omp parallel for
		for (size_t i = 0; i < primitive_count; ++i)
			primitives_copy[i] = primitives[indices[i]];
		return primitives_copy;
	}

	/// Compute bounding cylinders and centers
	template <typename Primitive, typename Scalar = typename Primitive::ScalarType>
	std::pair<std::unique_ptr<BoundingCyl<Scalar>[]>, std::unique_ptr<Vector3<Scalar>[]>>
		compute_bounding_cylinders_and_centers(const Primitive* primitives, size_t primitive_count)
	{
		auto bounding_boxes = std::make_unique<BoundingCyl<Scalar>[]>(primitive_count);
		auto centers = std::make_unique<Vector3<Scalar>[]>(primitive_count);

#pragma omp parallel for
		for (size_t i = 0; i < primitive_count; ++i) {
			bounding_boxes[i] = primitives[i].bounding_cyl();
			centers[i] = primitives[i].center();
		}

		return std::make_pair(std::move(bounding_boxes), std::move(centers));
	}

	/// Computes the bounding box and the center of each primitive in given array.
	template <typename Primitive, typename Scalar = typename Primitive::ScalarType>
	std::pair<std::unique_ptr<BoundingBox<Scalar>[]>, std::unique_ptr<Vector3<Scalar>[]>>
		compute_bounding_boxes_and_centers(const Primitive* primitives, size_t primitive_count)
	{
		auto bounding_boxes = std::make_unique<BoundingBox<Scalar>[]>(primitive_count);
		auto centers = std::make_unique<Vector3<Scalar>[]>(primitive_count);

#pragma omp parallel for
		for (size_t i = 0; i < primitive_count; ++i) {
			bounding_boxes[i] = primitives[i].bounding_box();
			centers[i] = primitives[i].center();
		}

		return std::make_pair(std::move(bounding_boxes), std::move(centers));
	}

	template <typename Scalar>
	BoundingCyl<Scalar> compute_bounding_cylinders_union(const BoundingCyl<Scalar>* bboxes, size_t count) {
		auto bbox = BoundingCyl<Scalar>::empty();

#pragma omp declare reduction \
        (bbox_extend:BoundingCyl<Scalar>:omp_out.extend(omp_in)) \
        initializer(omp_priv = BoundingCyl<Scalar>::empty())

#pragma omp parallel for reduction(bbox_extend: bbox)
		for (size_t i = 0; i < count; ++i)
			bbox.extend(bboxes[i]);

		std::cout << "global BB surface " << bbox.surface() << std::endl;
		return bbox;
	}

	/// Computes the union AABB of all the bounding cylinders in the given array.
	template <typename Scalar>
	BoundingBox<Scalar> compute_bounding_boxes_union(const BoundingCyl<Scalar>* bboxes, size_t count) {
		auto bbox = BoundingBox<Scalar>::empty();

#pragma omp declare reduction \
        (bbox_extend:BoundingBox<Scalar>:omp_out.extend(omp_in)) \
        initializer(omp_priv = BoundingBox<Scalar>::empty())

#pragma omp parallel for reduction(bbox_extend: bbox)
		for (size_t i = 0; i < count; ++i)
			bbox.extend(bboxes[i].AABB());

		std::cout << "global BB surface " << bbox.surface() << std::endl;
		return bbox;
	}

	/// Computes the union of all the bounding boxes in the given array.
	template <typename Scalar>
	BoundingBox<Scalar> compute_bounding_boxes_union(const BoundingBox<Scalar>* bboxes, size_t count) {
		auto bbox = BoundingBox<Scalar>::empty();

#pragma omp declare reduction \
        (bbox_extend:BoundingBox<Scalar>:omp_out.extend(omp_in)) \
        initializer(omp_priv = BoundingBox<Scalar>::empty())

#pragma omp parallel for reduction(bbox_extend: bbox)
		for (size_t i = 0; i < count; ++i)
			bbox.extend(bboxes[i]);

		std::cout << "global BB surface " << bbox.surface() << std::endl;
		return bbox;
	}

	/// Perpendicular distance of a point to a line
	template <typename Scalar>
	Scalar point2line(/*point:*/ Vector3<Scalar> p, /*point on line:*/ Vector3<Scalar> a, /*vector*/ Vector3<Scalar> v) {
		return length((p - a) - dot(p - a, v) * v);
	}

	/// Perpendicular projection of a point to a line
	template <typename Scalar>
	Vector3<Scalar> pointOnline(Vector3<Scalar> p, Vector3<Scalar> a, Vector3<Scalar> v) {
		return a + dot(p - a, v) * v;
	}

	/// 2 in 1 for the above, return value is the point, the distance is in the parameter dist
	template <typename Scalar>
	Vector3<Scalar> pointOn2line(Vector3<Scalar> p, Vector3<Scalar> a, Vector3<Scalar> v, Scalar& dist) {
		auto dotvv = dot(v, v);
		Vector3<Scalar> h = dotvv == 0 ? Scalar(0) * v : (dot(p - a, v) / dotvv) * v;
		dist = length((p - a) - h);
		return a + h;
	}

	/// Projects a cylinder onto a line, returning two endpoints of the projection, plus a radius
	template <typename Scalar>
	void projectCylOnLine(BoundingCyl<Scalar> cyl, Vector3<Scalar> axis, Vector3<Scalar> axpoint, Vector3<Scalar>& p1, Vector3<Scalar>& p2, Scalar& rad) {
		Scalar cosA, sinA, x, d1, d2;
		Vector3<Scalar> ax;
		// just in case axis was not normalized
		ax = normalize(axis);
		cosA = dot(ax, cyl.axis);
		sinA = 1 - cosA * cosA;
		sinA = sinA < 0 ? 0 : sqrt(sinA);
		x = sinA * cyl.r;

		p1 = pointOn2line(cyl.c, axpoint, ax, d1);
		p1 = p1 - ax * x;
		//if (isnan(p1[0]) != 0)
		//	p1 = Vector3<Scalar>(0);

		p2 = pointOn2line(cyl.c + cyl.axis * cyl.h, axpoint, ax, d2);
		p2 = p2 + ax * x;
		//if (isnan(p2[0]) != 0)
		//	p2 = Vector3<Scalar>(0);

		//if (isnan(d1) != 0)
		//	d1 = 0;
		//if (isnan(d2) != 0)
		//	d2 = 0;
		rad = cyl.r * cosA + std::max(d1, d2);
	}

	/// Establish the height and center of the constructed cylinder
	template <typename Scalar>
	Scalar getHeight(Vector3<Scalar> c1, Vector3<Scalar> c1h, Vector3<Scalar> c2, Vector3<Scalar> c2h, Vector3<Scalar>& center) {
		//Vector3<Scalar> points[4];
		//points[0] = c1;
		//points[1] = c1h;
		//points[2] = c2;
		//points[3] = c2h;
		//Scalar height = 0;
		//Scalar h;
		//for (size_t i = 0; i < 4; i++) {
		//	for (size_t j = 0; j < 4; j++) {
		//		if (i != j) {
		//			h = length(points[i] - points[j]);
		//			if (h > height) {
		//				height = h;
		//				center = points[i];
		//			}
		//		}
		//	}
		//}
		//return height;

		Scalar height = length(c1 - c1h);
		center = c1;
		Scalar hh = length(c1 - c2h);
		height = (hh > height) ? hh : height;
		hh = length(c2 - c1h);
		if (hh > height) {
			height = hh;
			center = c2;
		}
		hh = length(c2 - c2h);
		if (hh > height) {
			height = hh;
			center = c2;
		}
		return height;
	}

} // namespace bvh

#endif
