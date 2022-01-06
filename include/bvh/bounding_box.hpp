#ifndef BVH_BOUNDING_BOX_HPP
#define BVH_BOUNDING_BOX_HPP

#include "bvh/vector.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <cassert>
#include <iostream>
#include <chrono>

std::chrono::steady_clock myclock;

namespace bvh {

	/// A bounding box, represented with two extreme points.
	template <typename Scalar>
	struct BoundingBox {
		Vector3<Scalar> min, max;

		BoundingBox() = default;
		bvh__always_inline__ BoundingBox(const Vector3<Scalar>& v) : min(v), max(v) {}
		bvh__always_inline__ BoundingBox(const Vector3<Scalar>& min, const Vector3<Scalar>& max) : min(min), max(max) {}

		bvh__always_inline__ BoundingBox& shrink(const BoundingBox& bbox) {
			min = bvh::max(min, bbox.min);
			max = bvh::min(max, bbox.max);
			return *this;
		}

		bvh__always_inline__ BoundingBox& extend(const BoundingBox& bbox) {
			min = bvh::min(min, bbox.min);
			max = bvh::max(max, bbox.max);
			return *this;
		}

		bvh__always_inline__ BoundingBox& extend(const Vector3<Scalar>& v) {
			min = bvh::min(min, v);
			max = bvh::max(max, v);
			return *this;
		}

		bvh__always_inline__ Vector3<Scalar> diagonal() const {
			return max - min;
		}

		bvh__always_inline__ Vector3<Scalar> center() const {
			return (max + min) * Scalar(0.5);
		}

		bvh__always_inline__ Scalar half_area() const {
			auto d = diagonal();
			return (d[0] + d[1]) * d[2] + d[0] * d[1];
		}

		bvh__always_inline__ Scalar surface() const {
			return Scalar(2) * half_area();
		}

		bvh__always_inline__ Scalar volume() const {
			auto d = diagonal();
			return d[0] * d[1] * d[2];
		}

		bvh__always_inline__ size_t largest_axis() const {
			auto d = diagonal();
			size_t axis = 0;
			if (d[0] < d[1]) axis = 1;
			if (d[axis] < d[2]) axis = 2;
			return axis;
		}

		bvh__always_inline__ Scalar largest_extent() const {
			return diagonal()[largest_axis()];
		}

		bvh__always_inline__ bool is_contained_in(const BoundingBox& other) const {
			return
				max[0] <= other.max[0] && min[0] >= other.min[0] &&
				max[1] <= other.max[1] && min[1] >= other.min[1] &&
				max[2] <= other.max[2] && min[2] >= other.min[2];
		}

		bvh__always_inline__ static BoundingBox full() {
			return BoundingBox(
				Vector3<Scalar>(-std::numeric_limits<Scalar>::max()),
				Vector3<Scalar>(std::numeric_limits<Scalar>::max()));
		}

		bvh__always_inline__ static BoundingBox empty() {
			return BoundingBox(
				Vector3<Scalar>(std::numeric_limits<Scalar>::max()),
				Vector3<Scalar>(-std::numeric_limits<Scalar>::max()));
		}
	};

	/// Special bounding box - a cylinder.
	template <typename Scalar>
	struct BoundingCyl {
		Vector3<Scalar> c, axis;
		Scalar h, r;

		BoundingCyl() = default;
		bvh__always_inline__ BoundingCyl(const Vector3<Scalar> center, const Scalar height, const Scalar radius) :
			c(center), h(height), r(radius) {}
		bvh__always_inline__ BoundingCyl(const Vector3<Scalar> center, const Vector3<Scalar> ax, const Scalar height, const Scalar radius) :
			c(center), axis(ax), h(height), r(radius) {}

		void printme() const {
			std::cout << "center " << c[0] << " " << c[1] << " " << c[2] << std::endl;
			std::cout << "axis " << axis[0] << " " << axis[1] << " " << axis[2] << std::endl;
			std::cout << "height " << h << std::endl;
			std::cout << "radius " << r << std::endl << std::endl;
		}

		bvh__always_inline__ BoundingBox<Scalar> AABB() const {
			//Scalar xmin, xmax, ymin, ymax, zmin, zmax;
			Vector3<Scalar> upper(c + axis * h);
			Vector3<Scalar> mins(min(c, upper)), maxes(max(c, upper));

			Scalar kx, ky, kz, x2, y2, z2;
			x2 = c[0] - upper[0]; x2 *= x2;
			y2 = c[1] - upper[1]; y2 *= y2;
			z2 = c[2] - upper[2]; z2 *= z2;
			Scalar sum = x2 + y2 + z2;
			kx = sum == 0 ? 0 : Scalar(sqrt((y2 + z2) / sum));
			ky = sum == 0 ? 0 : Scalar(sqrt((x2 + z2) / sum));
			kz = sum == 0 ? 0 : Scalar(sqrt((x2 + y2) / sum));

			mins -= Vector3<Scalar>(kx * r, ky * r, kz * r);
			maxes += Vector3<Scalar>(kx * r, ky * r, kz * r);

			return BoundingBox(mins, maxes);
		}

		/// Make an AABB of the cylinder and use the min of that box
		bvh__always_inline__ Vector3<Scalar> getminpoint() const {
			return AABB().min;
		}

		bvh__always_inline__ Vector3<Scalar> upcenter() const {
			return c + normalize(axis) * h;
		}

		bvh__always_inline__ Vector3<Scalar> center() const {
			return c + (h * Scalar(0.5)) * axis;
		}

		bvh__always_inline__ Scalar surface() const {
			return Scalar(2) * M_PI * r * (h + r);
		}

		bvh__always_inline__ Scalar half_area() const {
			return Scalar(M_PI) * r * (h + r);
		}

		/// AABB diagonal
		bvh__always_inline__ Vector3<Scalar> diagAABB() const {
			return AABB().diagonal();
		}

		bvh__always_inline__ Vector3<Scalar> diagonal() const {
			Vector3<Scalar> radvec(axis[2], Scalar(0), -axis[0]);
			radvec = normalize(radvec) * r;
			Vector3<Scalar> p = (c + axis * h) + radvec;
			return p - (c - radvec);
		}

		bvh__always_inline__ Scalar volume() const {
			return M_PI * M_PI * r * h;
		}

		bvh__always_inline__ bool is_point_inside(Vector3<Scalar> p) {
			// 1. is its projection closer to the axis than R?
			Scalar d;
			Vector3<Scalar> proj = pointOn2line(p, this->c, this->axis, d);
			// if not, return false
			if (d > this->r)
				return false;

			// 2. if yes, is the point's projection between the two centers?
			// if not, return  false
			if (dot(this->axis, normalize(p - this->c)) != 1
				|| dot(this->axis, normalize(upcenter() - p)) != 1)
				return false;

			return true;
		}

		bvh__always_inline__ bool is_in_eps_range(Scalar a, Scalar b, Scalar e) {
			if (std::abs(a - b) > e) return false;
			return true;
		}

		bvh__always_inline__ bool is_around(Vector3<Scalar> a, Vector3<Scalar> b, Scalar e) {
			if (std::abs(a[0] - b[0]) > e) return false;
			if (std::abs(a[1] - b[1]) > e) return false;
			if (std::abs(a[2] - b[2]) > e) return false;
			return true;
		}


		bvh__always_inline__ bool is_inside(const BoundingCyl& other, BoundingCyl& big) {
			// 1. just project the smaller cylinder onto the bigger axis
			Vector3<Scalar> p1, p2;
			Scalar rad;
			BoundingCyl smaller = *this, bigger = other;
			if (other.volume() < this->volume()) {
				smaller = other;
				bigger = *this;
			}
			Scalar e = 0.00001; // epsilon

			projectCylOnLine(smaller, bigger.axis, bigger.c, p1, p2, rad);
			// 2. check if rad is smaller than the bigger radius
			if (rad > bigger.r)
				if (!is_in_eps_range(rad, bigger.r, e))
					return false;

			// 3. check if the two returned points are between the bigger axis endpoints
			// p1 - a, p1 - b
			if (dot(p1 - bigger.c, p1 - bigger.c + bigger.h * bigger.axis) > 0)
				if (!is_around(p1, bigger.c, e) && !is_around(p1, bigger.c + bigger.h * bigger.axis, e))
					return false;
			if (dot(p2 - bigger.c, p2 - bigger.c + bigger.h * bigger.axis) > 0)
				if (!is_around(p2, bigger.c, e) && !is_around(p2, bigger.c + bigger.h * bigger.axis, e))
					return false;

			// if all holds
			big = bigger;
			return true;
		}

		// Construction methods
		// "extend" the cylinder by another cylinder (enclose two cylinders by a new one)
		// basically compute the union cylinder of the existing cylinder and the new one
		bvh__always_inline__ BoundingCyl& extend(const BoundingCyl& bbox) {
			Vector3<Scalar> center(0);
			Scalar r1 = 0, r2 = 0, d1, d2, height = 0;
			Vector3<Scalar> A, B, C, D;

			// general position
			A = bbox.c;
			B = this->c;
			C = bbox.c + bbox.h * bbox.axis;
			D = this->c + this->h * this->axis;

			// vanishing cylinder problem
			// might be that their axes point towards each other (more than 120 deg.)
			// swap the endpoints and invert the axis direction
			if (dot(bbox.axis, this->axis) < -0.8) {
				std::swap(B, D);
				this->axis = normalize(D - B);
				this->c = B;
			}

			// test if it is already fully contained in the other
			BoundingCyl container;
			if (this->is_inside(bbox, container)) {
				return *this = container;
			}
			//----------------------------------------------------
			Vector3<Scalar> u, v, ax, mid1, mid2;
			// Use the weighted variant
			Scalar vol1 = bbox.volume();
			Scalar vol2 = this->volume();
			Scalar w = (vol1 + vol2 == 0) ? 0 : vol2 / (vol1 + vol2);
			w = w == 0 ? 0.5 : w;
			u = normalize(B - A);
			v = normalize(D - C);
			d1 = length(B - A) * w;
			d2 = length(D - C) * w;

			mid1 = A + d1 * u;
			mid2 = C + d2 * v;
			ax = normalize(mid2 - mid1);

			height = length(mid2 - mid1);
			Vector3<Scalar> c1, c1h, c2, c2h;

			projectCylOnLine(bbox, ax, mid1, c1, c1h, r1);
			projectCylOnLine(*this, ax, mid1, c2, c2h, r2);

			height = getHeight(c1, c1h, c2, c2h, center);
			BoundingCyl result;
			result.h = height;
			result.r = std::max(r1, r2);
			result.c = center;
			result.axis = ax;
			
			bool optim = false;
			// optimize - brute force find best from the 100 positions
			if (optim) {
				Scalar lenBA = length(B - A);
				Scalar lenDC = length(D - C);
				Scalar numsteps = 10;

				BoundingCyl helper;
				//Scalar height;
				for (int k = 0; k < numsteps; k++) { // numsteps steps on upper line
					mid1 = A + (lenBA * (Scalar)k / numsteps) * (u);
					for (int l = 0; l < numsteps; l++) { // numsteps steps on lower line
						mid2 = C + (lenDC * (Scalar)l / numsteps) * (v);

						ax = normalize(mid2 - mid1);
						projectCylOnLine(bbox, ax, mid1, c1, c1h, r1);
						projectCylOnLine(*this, ax, mid1, c2, c2h, r2);
						// compute height
						height = getHeight(c1, c1h, c2, c2h, center);

						helper.h = height;
						helper.r = std::max(r1, r2);
						helper.c = center;
						helper.axis = ax;

						if (helper.surface() < result.surface()) {
							result = helper;
						}
					}
				}
			}

			*this = result;
			return *this;
		}

		bvh__always_inline__ static BoundingCyl empty() {
			return BoundingCyl(
				Vector3<Scalar>(0),
				Vector3<Scalar>(0),
				0, 0);
		}
	};
} // namespace bvh

#endif
