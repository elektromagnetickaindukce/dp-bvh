#ifndef BVH_NODE_INTERSECTORS_HPP
#define BVH_NODE_INTERSECTORS_HPP

#include "bvh/vector.hpp"
#include "bvh/ray.hpp"
#include "bvh/platform.hpp"
#include "bvh/utilities.hpp"

namespace bvh {

	/// New node intersector method for cylinder shaped bounding boxes.
	template <typename Bvh>
	struct CustomNodeIntersector {
		using Scalar = typename Bvh::ScalarType;
		Vector3<Scalar> inverse_direction;
		Scalar invdiridx;

		CustomNodeIntersector(const Ray<Scalar>& ray) {
			inverse_direction = ray.direction.inverse();
			invdiridx = inverse_direction[0] == 0 ? (inverse_direction[1] == 0 ? (inverse_direction[2] == 0 ? -1 : 2) : 1) : 0;
		}

		bvh__always_inline__
		std::pair<Scalar, Scalar> intersect(const typename Bvh::CustomNode& node, const Ray<Scalar>& ray) const {
			// ray has origin and direction members
			// CustomNode holds a bounding cylinder

			// 1. step 
			// compute A, B, C
			Vector3<Scalar> d_p = ray.origin - node.p1;
			Scalar dot_vva = dot(node.axis, ray.direction);
			Scalar dot_dpva = dot(node.axis, d_p);
			Scalar r_2 = node.r * node.r;
			Scalar A, B, C;
			Vector3<Scalar> v, v2;
			v = ray.direction - dot_vva * node.axis;
			A = dot(v, v);
			v2 = d_p - dot_dpva * node.axis;
			B = Scalar(2) * dot(v, v2);
			C = dot(v2, v2) - r_2;

			// 2. step
			// first check
			Scalar sqrterm = B * B - Scalar(4) * A * C;
			if (sqrterm < 0)
				return std::pair<Scalar, Scalar>(
					std::numeric_limits<Scalar>::max(),
					-std::numeric_limits<Scalar>::max()); // no intersection

			Vector<Scalar, 4> tvec(std::numeric_limits<Scalar>::max());
			// 3. step
			// intersecting an infinite cylinder
			// solve for t1, t2
			Scalar t1 = (-B + sqrt(sqrterm)) / (Scalar(2) * A);
			Scalar t2 = (-B - sqrt(sqrterm)) / (Scalar(2) * A);
			Vector3<Scalar> p2 = node.p1 + node.h * node.axis;
			if (t1 > 0) {
				Vector3<Scalar> q1 = ray.origin + t1 * ray.direction;
				if (dot(node.axis, q1 - node.p1) > 0 &&
					dot(node.axis, q1 - p2) < 0)
					tvec[0] = t1;
			}
			if (t2 > 0) {
				Vector3<Scalar> q2 = ray.origin + t2 * ray.direction;
				if (dot(node.axis, q2 - node.p1) > 0 &&
					dot(node.axis, q2 - p2) < 0)
					tvec[1] = t2;
			}
			// 4. step
			// intersect with the cylinder caps
			// find t3, t4 if they exist
			if (invdiridx != -1) {
				Scalar t3 = (node.p1[invdiridx] - ray.origin[invdiridx]) * inverse_direction[invdiridx];
				Scalar t4 = (p2[invdiridx] - ray.origin[invdiridx]) * inverse_direction[invdiridx];
				Vector3<Scalar> q = ray.origin + t3 * ray.direction;
				Vector3<Scalar> qq = q - node.p1;
				if (t3 > 0 && dot(qq, qq) < r_2)
					tvec[2] = t3;
				q = ray.origin + t4 * ray.direction;
				qq = q - p2;
				if (t4 > 0 && dot(qq, qq) < r_2)
					tvec[3] = t4;
			}

			// 5. step
			// return smallest t from tvec
			return std::make_pair(
				robust_min(tvec[0], robust_min(tvec[1], robust_min(tvec[2], tvec[3]))),
				robust_max(tvec[0], robust_max(tvec[1], robust_max(tvec[2], tvec[3])))
			);
		}

	//protected:
		~CustomNodeIntersector() {}
	};

	/// Base class for ray-node intersection algorithms. Does ray octant classification.
	template <typename Bvh, typename Derived>
	struct NodeIntersector {
		using Scalar = typename Bvh::ScalarType;

		std::array<int, 3> octant;

		NodeIntersector(const Ray<Scalar>& ray)
			: octant{
				ray.direction[0] < Scalar(0),
				ray.direction[1] < Scalar(0),
				ray.direction[2] < Scalar(0)
		}
		{}

		bvh__always_inline__
			Scalar intersect_axis(int axis, const Vector3<Scalar>& p, const Ray<Scalar>& ray) const {
			return static_cast<const Derived*>(this)->intersect_axis(axis, p, ray);
		}

		bvh__always_inline__
			std::pair<Scalar, Scalar> intersect(const typename Bvh::Node& node, const Ray<Scalar>& ray) const {
			Vector3<Scalar> entry, exit;
			entry[0] = intersect_axis(0, node.bounds[0 * 2 + octant[0]], ray);
			entry[1] = intersect_axis(1, node.bounds[1 * 2 + octant[1]], ray);
			entry[2] = intersect_axis(2, node.bounds[2 * 2 + octant[2]], ray);
			exit[0] = intersect_axis(0, node.bounds[0 * 2 + 1 - octant[0]], ray);
			exit[1] = intersect_axis(1, node.bounds[1 * 2 + 1 - octant[1]], ray);
			exit[2] = intersect_axis(2, node.bounds[2 * 2 + 1 - octant[2]], ray);
			// Note: This order for the min/max operations is guaranteed not to produce NaNs
			return std::make_pair(
				robust_max(entry[0], robust_max(entry[1], robust_max(entry[2], ray.tmin))),
				robust_min(exit[0], robust_min(exit[1], robust_min(exit[2], ray.tmax))));
		}

	protected:
		~NodeIntersector() {}
	};

	/// Fully robust ray-node intersection algorithm (see "Robust BVH Ray Traversal", by T. Ize).
	template <typename Bvh>
	struct RobustNodeIntersector : public NodeIntersector<Bvh, RobustNodeIntersector<Bvh>> {
		using Scalar = typename Bvh::ScalarType;

		// Padded inverse direction to avoid false-negatives in the ray-node test.
		Vector3<Scalar> padded_inverse_direction;

		RobustNodeIntersector(const Ray<Scalar>& ray)
			: NodeIntersector<Bvh, RobustNodeIntersector<Bvh>>(ray)
		{
			auto inverse_direction = ray.direction.inverse();

			padded_inverse_direction = Vector3<Scalar>(
				add_ulp_magnitude(inverse_direction[0], 2),
				add_ulp_magnitude(inverse_direction[1], 2),
				add_ulp_magnitude(inverse_direction[2], 2));
		}

		bvh__always_inline__
			Scalar intersect_axis(int axis, const Vector3<Scalar>& p, const Ray<Scalar>& ray) const {
			return (p[axis] - ray.origin[axis]) * padded_inverse_direction[axis];
		}

		using NodeIntersector<Bvh, RobustNodeIntersector<Bvh>>::intersect;
	};

	/// Semi-robust, fast ray-node intersection algorithm.
	template <typename Bvh>
	struct FastNodeIntersector : public NodeIntersector<Bvh, FastNodeIntersector<Bvh>> {
		using Scalar = typename Bvh::ScalarType;

		Vector3<Scalar> scaled_origin;
		Vector3<Scalar> inverse_direction;

		FastNodeIntersector(const Ray<Scalar>& ray)
			: NodeIntersector<Bvh, FastNodeIntersector<Bvh>>(ray)
		{
			inverse_direction = ray.direction.inverse();
			scaled_origin = -ray.origin * inverse_direction;
		}

		bvh__always_inline__
			Scalar intersect_axis(int axis, const Vector3<Scalar>& p, const Ray<Scalar>&) const {
			return fast_multiply_add(p[axis], inverse_direction[axis], scaled_origin[axis]);
		}

		using NodeIntersector<Bvh, FastNodeIntersector<Bvh>>::intersect;
	};

} // namespace bvh

#endif
