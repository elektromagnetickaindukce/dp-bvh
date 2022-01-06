#ifndef BVH_SINGLE_RAY_TRAVERSAL_HPP
#define BVH_SINGLE_RAY_TRAVERSAL_HPP

#include <cassert>

#include "bvh/bvh.hpp"
#include "bvh/ray.hpp"
#include "bvh/node_intersectors.hpp"
#include "bvh/utilities.hpp"

namespace bvh {

	/// Single ray traversal algorithm, using the provided ray-node intersector.
	template <typename Bvh, typename NodeIntersector = FastNodeIntersector<Bvh>, size_t StackSize = 256>
	class SingleRayTraverser {
	public:
		static constexpr size_t stack_size = StackSize;

	private:
		using Scalar = typename Bvh::ScalarType;

		struct Stack {
			using Element = const typename Bvh::Node*;

			Element elements[stack_size];
			size_t size = 0;

			void push(const Element& t) {
				assert(size < stack_size);
				elements[size++] = t;
			}

			Element pop() {
				assert(!empty());
				return elements[--size];
			}

			bool empty() const { return size == 0; }
		};

		struct StackC {
			using Element = const typename Bvh::CustomNode*;

			Element elements[stack_size];
			size_t size = 0;

			void push(const Element& t) {
				assert(size < stack_size);
				elements[size++] = t;
			}

			Element pop() {
				assert(!empty());
				return elements[--size];
			}

			bool empty() const { return size == 0; }
		};

		/// New leaf intersector variant
		template <typename PrimitiveIntersector, typename Statistics>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>& intersect_leaf(
				const typename Bvh::CustomNode& node,
				Ray<Scalar>& ray,
				std::optional<typename PrimitiveIntersector::Result>& best_hit,
				PrimitiveIntersector& primitive_intersector,
				Statistics& statistics) const
		{
			assert(node.is_leaf);
			size_t begin = node.first_child_or_primitive;
			size_t end = begin + node.primitive_count;
			statistics.intersections += end - begin;
			for (size_t i = begin; i < end; ++i) {
				if (auto hit = primitive_intersector.intersect(i, ray)) {
					best_hit = hit;
					if (primitive_intersector.any_hit)
						return best_hit;
					ray.tmax = hit->distance();
				}
			}
			return best_hit;
		}

		template <typename PrimitiveIntersector, typename Statistics>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>& intersect_leaf(
				const typename Bvh::Node& node,
				Ray<Scalar>& ray,
				std::optional<typename PrimitiveIntersector::Result>& best_hit,
				PrimitiveIntersector& primitive_intersector,
				Statistics& statistics) const
		{
			assert(node.is_leaf);
			size_t begin = node.first_child_or_primitive;
			size_t end = begin + node.primitive_count;
			statistics.intersections += end - begin;
			for (size_t i = begin; i < end; ++i) {
				if (auto hit = primitive_intersector.intersect(i, ray)) {
					best_hit = hit;
					if (primitive_intersector.any_hit)
						return best_hit;
					ray.tmax = hit->distance();
				}
			}
			return best_hit;
		}

		/// intersect function for cylinders 
		template <typename PrimitiveIntersector, typename Statistics>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			intersectC(Ray<Scalar> ray, PrimitiveIntersector& primitive_intersector, Statistics& statistics) const {
			auto best_hit = std::optional<typename PrimitiveIntersector::Result>(std::nullopt);

			// If the root is a leaf, intersect it and return
			if (bvh__unlikely(bvh.cnodes[0].is_leaf))
				return intersect_leaf(bvh.cnodes[0], ray, best_hit, primitive_intersector, statistics);

			bvh::CustomNodeIntersector<Bvh> node_intersector(ray);

			// This traversal loop is eager, because it immediately processes leaves instead of pushing them on the stack.
			// This is generally beneficial for performance because intersections will likely be found which will
			// allow to cull more subtrees with the ray-box test of the traversal loop.
			StackC stack;
			const auto* node = bvh.cnodes.get();
			while (true) {
				statistics.traversal_steps++;

				auto first_child = node->first_child_or_primitive;
				const auto* left_child = &bvh.cnodes[first_child + 0];
				const auto* right_child = &bvh.cnodes[first_child + 1];
				auto distance_left = node_intersector.intersect(*left_child, ray);
				auto distance_right = node_intersector.intersect(*right_child, ray);

				if (distance_left.first <= distance_left.second) {
					if (bvh__unlikely(left_child->is_leaf)) {
						if (intersect_leaf(*left_child, ray, best_hit, primitive_intersector, statistics) &&
							primitive_intersector.any_hit)
							break;
						left_child = nullptr;
					}
				}
				else
					left_child = nullptr;

				if (distance_right.first <= distance_right.second) {
					if (bvh__unlikely(right_child->is_leaf)) {
						if (intersect_leaf(*right_child, ray, best_hit, primitive_intersector, statistics) &&
							primitive_intersector.any_hit)
							break;
						right_child = nullptr;
					}
				}
				else
					right_child = nullptr;

				if (bvh__likely((left_child != NULL) ^ (right_child != NULL))) {
					node = left_child != NULL ? left_child : right_child;
				}
				else if (bvh__unlikely((left_child != NULL) & (right_child != NULL))) {
					if (distance_left.first > distance_right.first)
						std::swap(left_child, right_child);
					stack.push(right_child);
					node = left_child;
				}
				else {
					if (stack.empty())
						break;
					node = stack.pop();
				}
			}

			return best_hit;
		}

		template <typename PrimitiveIntersector, typename Statistics>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			intersectH(Ray<Scalar> ray, PrimitiveIntersector& primitive_intersector, Statistics& statistics) const {
			auto best_hit = std::optional<typename PrimitiveIntersector::Result>(std::nullopt);

			// If the root is a leaf, intersect it and return
			if (bvh__unlikely(bvh.nodes[0].is_leaf))
				return intersect_leaf(bvh.nodes[0], ray, best_hit, primitive_intersector, statistics);

			NodeIntersector node_intersector(ray);
			bvh::CustomNodeIntersector<Bvh> cnode_intersector(ray);

			Stack stack;
			StackC stackc; // second stack for cylinder nodes
			const auto* node = bvh.nodes.get();
			const auto* cnode = bvh.cnodes.get();
			while (true) {
				statistics.traversal_steps++;

				if (node->is_leaf) {
					cnode = &bvh.cnodes[node->origin];

					// if cylinder is a leaf
					if (cnode->is_leaf) {
						if (intersect_leaf(*cnode, ray, best_hit, primitive_intersector, statistics) &&
							primitive_intersector.any_hit)
							if (stack.empty())
								break;
						node = stack.pop();
						continue;
					}

					while (true) {
						statistics.traversal_steps++;
						auto first_child = cnode->first_child_or_primitive;
						const auto* left = &bvh.cnodes[first_child + 0];
						const auto* right = &bvh.cnodes[first_child + 1];
						auto distance_left = cnode_intersector.intersect(*left, ray);
						auto distance_right = cnode_intersector.intersect(*right, ray);

						if (distance_left.first <= distance_left.second) {
							if (bvh__unlikely(left->is_leaf)) {
								if (intersect_leaf(*left, ray, best_hit, primitive_intersector, statistics) &&
									primitive_intersector.any_hit)
									break;
								left = nullptr;
							}
						}
						else
							left = nullptr;

						if (distance_right.first <= distance_right.second) {
							if (bvh__unlikely(right->is_leaf)) {
								if (intersect_leaf(*right, ray, best_hit, primitive_intersector, statistics) &&
									primitive_intersector.any_hit)
									break;
								right = nullptr;
							}
						}
						else
							right = nullptr;

						// push to the cylinder stack
						if (bvh__likely((left != NULL) ^ (right != NULL))) {
							cnode = left != NULL ? left : right;
						}
						else if (bvh__unlikely((left != NULL) & (right != NULL))) {
							if (distance_left.first > distance_right.first)
								std::swap(left, right);
							stackc.push(right);
							cnode = left;
						}
						else {
							if (stackc.empty())
								break;
							cnode = stackc.pop();
						}
					}
					if (stack.empty())
						break;
					node = stack.pop();
					continue;
				}

				auto first_child = node->first_child_or_primitive;
				const auto* left_child = &bvh.nodes[first_child + 0];
				const auto* right_child = &bvh.nodes[first_child + 1];
				auto distance_left = node_intersector.intersect(*left_child, ray);
				auto distance_right = node_intersector.intersect(*right_child, ray);

				if (distance_left.first <= distance_left.second) {
					// if left child is leaf, change to cylinders
					if (bvh__unlikely(left_child->is_leaf)) {
						const auto* cleaf = &bvh.cnodes[left_child->origin];
						if (cleaf->is_leaf) {
							if (intersect_leaf(*cleaf, ray, best_hit, primitive_intersector, statistics) &&
								primitive_intersector.any_hit)
								break;
							left_child = nullptr;
						}
					}
				}
				else
					left_child = nullptr;

				if (distance_right.first <= distance_right.second) {
					// if right child is leaf, change to cylinders
					if (bvh__unlikely(right_child->is_leaf)) {
						const auto* cleaf = &bvh.cnodes[right_child->origin];
						if (cleaf->is_leaf) {
							if (intersect_leaf(*cleaf, ray, best_hit, primitive_intersector, statistics) &&
								primitive_intersector.any_hit)
								break;
							right_child = nullptr;
						}
					}
				}
				else
					right_child = nullptr;

				if (bvh__likely((left_child != NULL) ^ (right_child != NULL))) {
					node = left_child != NULL ? left_child : right_child;
				}
				else if (bvh__unlikely((left_child != NULL) & (right_child != NULL))) {
					if (distance_left.first > distance_right.first)
						std::swap(left_child, right_child);
					stack.push(right_child);
					node = left_child;
				}
				else {
					if (stack.empty())
						break;
					node = stack.pop();
				}
			}

			return best_hit;
		}

		template <typename PrimitiveIntersector, typename Statistics>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			intersect(Ray<Scalar> ray, PrimitiveIntersector& primitive_intersector, Statistics& statistics) const {
			auto best_hit = std::optional<typename PrimitiveIntersector::Result>(std::nullopt);

			// If the root is a leaf, intersect it and return
			if (bvh__unlikely(bvh.nodes[0].is_leaf))
				return intersect_leaf(bvh.nodes[0], ray, best_hit, primitive_intersector, statistics);

			NodeIntersector node_intersector(ray);

			// This traversal loop is eager, because it immediately processes leaves instead of pushing them on the stack.
			// This is generally beneficial for performance because intersections will likely be found which will
			// allow to cull more subtrees with the ray-box test of the traversal loop.
			Stack stack;
			const auto* node = bvh.nodes.get();
			while (true) {
				statistics.traversal_steps++;

				auto first_child = node->first_child_or_primitive;
				const auto* left_child = &bvh.nodes[first_child + 0];
				const auto* right_child = &bvh.nodes[first_child + 1];
				auto distance_left = node_intersector.intersect(*left_child, ray);
				auto distance_right = node_intersector.intersect(*right_child, ray);

				if (distance_left.first <= distance_left.second) {
					if (bvh__unlikely(left_child->is_leaf)) {
						if (intersect_leaf(*left_child, ray, best_hit, primitive_intersector, statistics) &&
							primitive_intersector.any_hit)
							break;
						left_child = nullptr;
					}
				}
				else
					left_child = nullptr;

				if (distance_right.first <= distance_right.second) {
					if (bvh__unlikely(right_child->is_leaf)) {
						if (intersect_leaf(*right_child, ray, best_hit, primitive_intersector, statistics) &&
							primitive_intersector.any_hit)
							break;
						right_child = nullptr;
					}
				}
				else
					right_child = nullptr;

				if (bvh__likely((left_child != NULL) ^ (right_child != NULL))) {
					node = left_child != NULL ? left_child : right_child;
				}
				else if (bvh__unlikely((left_child != NULL) & (right_child != NULL))) {
					if (distance_left.first > distance_right.first)
						std::swap(left_child, right_child);
					stack.push(right_child);
					node = left_child;
				}
				else {
					if (stack.empty())
						break;
					node = stack.pop();
				}
			}

			return best_hit;
		}

		const Bvh& bvh;

	public:
		/// Statistics collected during traversal.
		struct Statistics {
			size_t traversal_steps = 0;
			size_t intersections = 0;
		};

		SingleRayTraverser(const Bvh& bvh)
			: bvh(bvh)
		{}

		/// Intersects the BVH with the given ray and intersector.
		template <typename PrimitiveIntersector>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			traverse(const Ray<Scalar>& ray, PrimitiveIntersector& intersector) const {
			struct {
				struct Empty {
					Empty& operator ++ (int) { return *this; }
					Empty& operator ++ () { return *this; }
					Empty& operator += (size_t) { return *this; }
				} traversal_steps, intersections;
			} statistics;
			return intersect(ray, intersector, statistics);
		}

		template <typename PrimitiveIntersector>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			traverse(const Ray<Scalar>& ray, PrimitiveIntersector& intersector, bool cyl, bool hybrid) const {
			struct {
				struct Empty {
					Empty& operator ++ (int) { return *this; }
					Empty& operator ++ () { return *this; }
					Empty& operator += (size_t) { return *this; }
				} traversal_steps, intersections;
			} statistics;
			return hybrid ? intersectH(ray, intersector, statistics) : cyl ? intersectC(ray, intersector, statistics) : intersect(ray, intersector, statistics);
		}

		/// Intersects the BVH with the given ray and intersector.
		/// Record statistics on the number of traversal and intersection steps.
		template <typename PrimitiveIntersector>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			traverse(const Ray<Scalar>& ray, PrimitiveIntersector& primitive_intersector, Statistics& statistics) const {
			return intersect(ray, primitive_intersector, statistics);
		}

		template <typename PrimitiveIntersector>
		bvh__always_inline__
			std::optional<typename PrimitiveIntersector::Result>
			traverse(const Ray<Scalar>& ray, PrimitiveIntersector& primitive_intersector, bool cyl, bool hybrid, Statistics& statistics) const {
			return hybrid ? intersectH(ray, primitive_intersector, statistics) : cyl ? intersectC(ray, primitive_intersector, statistics) : intersect(ray, primitive_intersector, statistics);
		}
	};

} // namespace bvh

#endif
