#ifndef BVH_TRIANGLE_HPP
#define BVH_TRIANGLE_HPP

#include <optional>
#include <cassert>
#include <math.h>
#include <limits>

#include "bvh/utilities.hpp"
#include "bvh/vector.hpp"
#include "bvh/bounding_box.hpp"
#include "bvh/ray.hpp"

namespace bvh {

	/// Triangle primitive, defined by three points, and using the Moeller-Trumbore test.
	/// By default, the normal is left-handed, which minimizes the number of operations in
	/// the intersection routine.
	template <typename Scalar, bool LeftHandedNormal = true>
	struct Triangle {
		struct Intersection {
			Scalar t, u, v;
			Scalar distance() const { return t; }
		};

		using ScalarType = Scalar;
		using IntersectionType = Intersection;

		Vector3<Scalar> p0, e1, e2, n;

		Triangle() = default;
		Triangle(const Vector3<Scalar>& p0, const Vector3<Scalar>& p1, const Vector3<Scalar>& p2)
			: p0(p0), e1(p0 - p1), e2(p2 - p0)
		{
			n = LeftHandedNormal ? cross(e1, e2) : cross(e2, e1);
		}

		Vector3<Scalar> p1() const { return p0 - e1; }
		Vector3<Scalar> p2() const { return p0 + e2; }

		/// Compute bounding cylinder
		BoundingCyl<Scalar> bounding_cyl() const {
			BoundingCyl<Scalar> box;
			Vector3<Scalar> A, B, C, mid;
			Vector3<Scalar> v = p1() - p2();
			Scalar l = 0;
			if (length(e1) <= length(e2)) {  // find shortest side
				A = p1(); B = p0; C = p2(); l = length(e1);
			}
			else {
				A = p0; B = p2(); C = p1(); l = length(e2);
			}
			if (length(v) <= l) {
				A = p2(); B = p1(); C = p0; l = length(v);
			}
			v = normalize(B - A);
			//assert(l != 0);
			//if (length(v) == 0) {
			//	v = Vector3<Scalar>(0);
			//	l = 0;
			//}
			//else {
			//	v = normalize(v);
			//}

			//if (isnan(v[0]) != 0 || isnan(v[1]) != 0 || isnan(v[2]) != 0) {
			//	v = Vector3<Scalar>(0); //(std::numeric_limits<Scalar>::min());
			//	l = 0; // length(v);
			//}
			//Scalar myeps = 0.0001f;
			//if (l < myeps)
			//	l = myeps;
			mid = A + (l / Scalar(2)) * v;
			Scalar l2 = 0;
			Vector3<Scalar> c1 = pointOn2line(A, C, mid - C, l);
			Vector3<Scalar> c2 = pointOn2line(B, C, mid - C, l);
			// reuse A and B
			//A = C - c1; B = C - c2;
			A = c1 - C; B = c2 - C;
			if (length(A) >= length(B)) {
				box.axis = normalize(A);
				//box.c = c1;
				box.h = length(A);
				//box.r = l;
			}
			else {
				box.axis = normalize(B);
				//box.c = c2;
				box.h = length(B);
				//box.r = l2;
			}
			box.r = l; //std::max(l, l2);
			box.c = C;
			//if (box.h < myeps)
			//	box.h = myeps;
			//if (box.r < myeps)
			//	box.r = myeps;

			// check for NaNs
			assert(isnan(box.r) == 0);
			//assert(isnan(box.h) == 0);

			//assert(isnan(box.c[0]) == 0);
			//assert(isnan(box.c[1]) == 0);
			//assert(isnan(box.c[2]) == 0);

			//assert(isnan(box.axis[0]) == 0);
			//assert(isnan(box.axis[1]) == 0);
			//assert(isnan(box.axis[2]) == 0);

			return box;
		}

		BoundingBox<Scalar> bounding_box() const {
			BoundingBox<Scalar> bbox(p0);
			bbox.extend(p1());
			bbox.extend(p2());
			return bbox;
		}

		Vector3<Scalar> center() const {
			return (p0 + p1() + p2()) * (Scalar(1.0) / Scalar(3.0));
		}

		std::pair<Vector3<Scalar>, Vector3<Scalar>> edge(size_t i) const {
			assert(i < 3);
			Vector3<Scalar> p[] = { p0, p1(), p2() };
			return std::make_pair(p[i], p[(i + 1) % 3]);
		}

		Scalar area() const {
			return length(n) * Scalar(0.5);
		}

		std::pair<BoundingBox<Scalar>, BoundingBox<Scalar>> split(size_t axis, Scalar position) const {
			Vector3<Scalar> p[] = { p0, p1(), p2() };
			auto left = BoundingBox<Scalar>::empty();
			auto right = BoundingBox<Scalar>::empty();
			auto split_edge = [=](const Vector3<Scalar>& a, const Vector3<Scalar>& b) {
				auto t = (position - a[axis]) / (b[axis] - a[axis]);
				return a + t * (b - a);
			};
			auto q0 = p[0][axis] <= position;
			auto q1 = p[1][axis] <= position;
			auto q2 = p[2][axis] <= position;
			if (q0) left.extend(p[0]);
			else    right.extend(p[0]);
			if (q1) left.extend(p[1]);
			else    right.extend(p[1]);
			if (q2) left.extend(p[2]);
			else    right.extend(p[2]);
			if (q0 ^ q1) {
				auto m = split_edge(p[0], p[1]);
				left.extend(m);
				right.extend(m);
			}
			if (q1 ^ q2) {
				auto m = split_edge(p[1], p[2]);
				left.extend(m);
				right.extend(m);
			}
			if (q2 ^ q0) {
				auto m = split_edge(p[2], p[0]);
				left.extend(m);
				right.extend(m);
			}
			return std::make_pair(left, right);
		}

		std::optional<Intersection> intersect(const Ray<Scalar>& ray) const {
			auto negate_when_right_handed = [](Scalar x) { return LeftHandedNormal ? x : -x; };

			auto c = p0 - ray.origin;
			auto r = cross(ray.direction, c);
			auto inv_det = negate_when_right_handed(1.0) / dot(n, ray.direction);

			auto u = dot(r, e2) * inv_det;
			auto v = dot(r, e1) * inv_det;
			auto w = Scalar(1.0) - u - v;

			// These comparisons are designed to return false
			// when one of t, u, or v is a NaN
			if (u >= 0 && v >= 0 && w >= 0) {
				auto t = negate_when_right_handed(dot(n, c)) * inv_det;
				if (t >= ray.tmin && t < ray.tmax)
					return std::make_optional(Intersection{ t, u, v });
			}

			return std::nullopt;
		}
	};

} // namespace bvh

#endif
