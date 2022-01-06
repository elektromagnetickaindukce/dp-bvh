#ifndef BVH_OBJ_EXPORTER_HPP
#define BVH_OBJ_EXPORTER_HPP

#include <queue>
#include <string>
#include <fstream>
#include <iostream>
#include "bvh/bvh.hpp"
#include "bvh/bounding_box.hpp"

//#define _USE_MATH_DEFINES
//#include <math.h>

namespace bvh {

	template <typename Bvh>
	class ObjExporter {
	public:
		using Scalar = typename Bvh::ScalarType;
		using Element = typename Bvh::CustomNode*;
		using Elem = typename Bvh::Node*;
		const Bvh& bvh;
		const std::string infile;
		const int numpoints;

	public:
		ObjExporter(const Bvh& bvh, int npoints = 20)
			: bvh(bvh), numpoints(npoints)
		{}

		ObjExporter(const Bvh& bvh, std::string filename, int npoints = 20)
			: bvh(bvh), numpoints(npoints), infile(filename.substr(0, filename.find_last_of(".")))
		{
			//infile = filename.substr(0, filename.find_last_of(".")); // remove the file extension
		}

		/// Export a bounding box
		void exportBox(typename Bvh::Node& box, std::ofstream& file, int id) {
			// bounds = min, max, min, max, min, max
			file << "v " << box.bounds[0] << " " << box.bounds[2] << " " << box.bounds[4] << std::endl;
			file << "v " << box.bounds[0] << " " << box.bounds[2] << " " << box.bounds[5] << std::endl;
			file << "v " << box.bounds[1] << " " << box.bounds[2] << " " << box.bounds[5] << std::endl;
			file << "v " << box.bounds[1] << " " << box.bounds[2] << " " << box.bounds[4] << std::endl;

			file << "v " << box.bounds[0] << " " << box.bounds[3] << " " << box.bounds[4] << std::endl;
			file << "v " << box.bounds[0] << " " << box.bounds[3] << " " << box.bounds[5] << std::endl;
			file << "v " << box.bounds[1] << " " << box.bounds[3] << " " << box.bounds[5] << std::endl;
			file << "v " << box.bounds[1] << " " << box.bounds[3] << " " << box.bounds[4] << std::endl << std::endl;

			file << "g cylinder" << std::to_string(id) << std::endl;
			file << "f " << "-8 -7 -6 -5\n";
			file << "f " << "-4 -3 -2 -1\n";

			file << "f " << "-1 -2 -6 -5\n";
			file << "f " << "-4 -3 -7 -8\n";

			file << "f " << "-2 -3 -7 -6\n";
			file << "f " << "-1 -4 -8 -5\n\n";
		}

		/// Export a cylinder
		//template<typename Bvh>
		void exportBox(typename Bvh::CustomNode& cyl, std::ofstream& file, int id, int npoints = 20) {
			// do the vertices first
			auto point = cyl.p1 + cyl.axis * cyl.h;
			file << "v " << cyl.p1[0] << " " << cyl.p1[1] << " " << cyl.p1[2] << std::endl; // 1
			file << "v " << point[0] << " " << point[1] << " " << point[2] << std::endl; // 2
			//file << std::endl;

			auto A = normalize(Vector3<Scalar>(cyl.axis[2], Scalar(0), -cyl.axis[0]));
			auto B = normalize(cross(A, cyl.axis));

			// add first two points
			auto pp = cyl.p1 + (A * Scalar(cos(0)) + B * Scalar(sin(0))) * cyl.r;
			file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
			pp = pp + cyl.axis * cyl.h;
			file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;

			for (int i = 0; i < npoints; i++) {
				auto v = Scalar(i) / npoints;
				pp = cyl.p1 + (A * Scalar(cos(2 * M_PI * v)) + B * Scalar(sin(2 * M_PI * v))) * cyl.r;
				file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
				pp = pp + cyl.axis * cyl.h;
				file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
			}

			file << std::endl;
			file << "g cylinder" << std::to_string((int)id) << std::endl;
			file << "usemtl cylinder" << std::endl;
			// then add the faces
			int i;
			int n = npoints * 2 + 4 + 1;
			for (i = 3; i <= npoints * 2 + 1; i += 2) {
				file << "f " << 1 - n << " " << i - n << " " << (i + 2) - n << std::endl;
				file << "f " << i - n << " " << (i + 1) - n << " " << (i + 2) - n << std::endl;
				file << "f " << (i + 2) - n << " " << (i + 1) - n << " " << (i + 3) - n << std::endl;
				file << "f " << (i + 3) - n << " " << (i + 1) - n << " " << 2 - n << std::endl;
			}
			// close the circle
			file << "f " << 1 - n << " " << i - n << " " << 3 - n << std::endl;
			file << "f " << i - n << " " << (i + 1) - n << " " << 3 - n << std::endl;
			file << "f " << 3 - n << " " << (i + 1) - n << " " << 4 - n << std::endl;
			file << "f " << 4 - n << " " << (i + 1) - n << " " << 2 - n << std::endl << std::endl;

		}

		void exportBox(typename BoundingCyl<Scalar> cyl, std::ofstream& file, int id, int npoints = 20) {
			// do the vertices first
			auto point = cyl.c + cyl.axis * cyl.h;
			file << "v " << cyl.c[0] << " " << cyl.c[1] << " " << cyl.c[2] << std::endl; // 1
			file << "v " << point[0] << " " << point[1] << " " << point[2] << std::endl; // 2

			auto A = normalize(Vector3<Scalar>(cyl.axis[2], Scalar(0), -cyl.axis[0]));
			auto B = normalize(cross(A, cyl.axis));

			// add first two points
			auto pp = cyl.c + (A * Scalar(cos(0)) + B * Scalar(sin(0))) * cyl.r;
			file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
			pp = pp + cyl.axis * cyl.h;
			file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;

			for (int i = 0; i < npoints; i++) {
				auto v = Scalar(i) / npoints;
				pp = cyl.c + (A * Scalar(cos(2 * M_PI * v)) + B * Scalar(sin(2 * M_PI * v))) * cyl.r;
				file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
				pp = pp + cyl.axis * cyl.h;
				file << "v " << pp[0] << " " << pp[1] << " " << pp[2] << std::endl;
			}

			file << std::endl;
			file << "g cylinder" << std::to_string((int)id) << std::endl;
			file << "usemtl cylinder" << std::endl;
			// then add the faces
			int i;
			int n = npoints * 2 + 4 + 1;
			for (i = 3; i <= npoints * 2 + 1; i += 2) {
				file << "f " << 1 - n << " " << i - n << " " << (i + 2) - n << std::endl;
				file << "f " << i - n << " " << (i + 1) - n << " " << (i + 2) - n << std::endl;
				file << "f " << (i + 2) - n << " " << (i + 1) - n << " " << (i + 3) - n << std::endl;
				file << "f " << (i + 3) - n << " " << (i + 1) - n << " " << 2 - n << std::endl;
			}
			// close the circle
			file << "f " << 1 - n << " " << i - n << " " << 3 - n << std::endl;
			file << "f " << i - n << " " << (i + 1) - n << " " << 3 - n << std::endl;
			file << "f " << 3 - n << " " << (i + 1) - n << " " << 4 - n << std::endl;
			file << "f " << 4 - n << " " << (i + 1) - n << " " << 2 - n << std::endl << std::endl;

		}

		/// Export CustomNode clusters
		void exportToFile(std::string infile, std::string filename, size_t level,
			std::unique_ptr<typename Bvh::CustomNode[]>& nodes, size_t begin, size_t end, Scalar& ret)
		{
			std::ofstream file;
			file.open(filename + "_" + std::to_string(level) + ".obj");
			file << infile << std::endl;
			std::ofstream specfile;
			int id = 0;
			Scalar surf = 0;
			for (size_t i = begin; i < end; ++i) {
				exportBox(nodes[i], file, id++);
				surf += nodes[i].bounding_box_proxy().to_bounding_box().surface();
			}
			ret = surf;
			file.close();
		}

		/// Export Node clusters
		void exportToFile(std::string infile, std::string filename, size_t level,
			std::unique_ptr<typename Bvh::Node[]>& nodes, size_t begin, size_t end, Scalar& ret)
		{
			std::ofstream file;
			file.open(filename + "_" + std::to_string(level) + ".obj");
			file << infile << std::endl;
			int id = 0;
			Scalar surf = 0;
			for (size_t i = begin; i < end; ++i) {
				exportBox(nodes[i], file, id++);
				surf += nodes[i].bounding_box_proxy().to_bounding_box().surface();
			}
			ret = surf;
			file.close();
		}

		void traverseExportHybrid() {
			///traverse the upper box part and then find the way to the cylinder part
			// q for the boxes
			std::queue<std::pair<Elem, int>> q;
			// second q for cylinders
			std::queue<std::pair<Element, int>> q2;

			int level = 0;
			q.push(std::make_pair(bvh.nodes.get(), level));
			auto node = q.front();

			// open files
			//std::ofstream bigstat;
			std::ofstream file;
			//bigstat.open(infile + ".txt", std::ios_base::out | std::ios_base::app);
			//file.open(infile + "_" + std::to_string(level) + ".obj");

			//file << str;
			//file << std::endl;

			int id = 0;

			// count the overall surface area on each level
			Scalar surface = 0;
			Scalar cumsum = 0;
			// and write into a separate file
			std::ofstream statfile;
			statfile.open("stats_box.txt");

			Scalar h;
			while (!q.empty()) {
				// pop the first
				node = q.front();
				q.pop();

				// process q2 until the level stays the same
				// or q2 is not empty
				while (!q2.empty()) {
					if (q2.front().second == level) {
						auto q2node = q2.front();
						q2.pop();
						//exportBox(*q2node.first, file, id++);
						h = (*q2node.first).bounding_box_proxy().to_bounding_box().surface();
						surface += h;
						cumsum += h;

						auto first_child = q2node.first->first_child_or_primitive;
						if (q2node.first->is_leaf) {
							//typename Bvh::IndexType i;
							//for (i = 0; i < q2node.first->primitive_count; i++) {
							//	exportBox(bvh.cnodes[first_child + i], file, id++);
							//}
						}
						else {
							auto* left_child = &bvh.cnodes[first_child + 0];
							auto* right_child = &bvh.cnodes[first_child + 1];

							if (left_child != NULL)
								q2.push(std::make_pair(left_child, q2node.second + 1));
							if (right_child != NULL)
								q2.push(std::make_pair(right_child, q2node.second + 1));
						}
					}
					else
						break;
				}

				// do something to the element
				// while the element level is equal to level, write to open file,
				// when it changes, close it and open a new file
				if (node.second != level) {
					// when changing the level, write the summed surface area into the statfile
					statfile << level << " " << (surface) << " " << cumsum << std::endl;
					surface = 0;

					file.close();
					level++;
					//file.open(infile + std::to_string((level)) + ".obj");

					//file << std::endl;
				}
				// export the node
				//exportBox(*node.first, file, id++);
				// add the surface 
				h = (*node.first).bounding_box_proxy().to_bounding_box().surface();
				surface += h;
				cumsum += h;

				// enqueue children, if they exist
				// if not leaf, add children to the queue
				auto first_child = node.first->first_child_or_primitive;

				if (node.first->is_leaf) {
					auto* cyl = &bvh.cnodes[node.first->origin];

					auto first_child_cyl = cyl->first_child_or_primitive;
					if (cyl->is_leaf) {
						//typename Bvh::IndexType i;
						//for (i = 0; i < cyl->primitive_count; i++) {
						//	exportBox(bvh.cnodes[first_child_cyl + i], file, id++);
						//}
					}
					else {
						auto* left_child = &bvh.cnodes[first_child_cyl + 0];
						auto* right_child = &bvh.cnodes[first_child_cyl + 1];

						if (left_child != NULL)
							q2.push(std::make_pair(left_child, node.second + 1));
						if (right_child != NULL)
							q2.push(std::make_pair(right_child, node.second + 1));
					}
				}
				else {
					auto* left_child = &bvh.nodes[first_child + 0];
					auto* right_child = &bvh.nodes[first_child + 1];

					if (left_child != NULL)
						q.push(std::make_pair(left_child, node.second + 1));
					if (right_child != NULL)
						q.push(std::make_pair(right_child, node.second + 1));
				}
			}
			// now process what's left in q2
			assert(!q2.empty());
			auto node2 = q2.front();

			// now empty the q2
			while (!q2.empty()) {
				node2 = q2.front();
				q2.pop();

				if (node2.second != level) {
					statfile << level << " " << surface << " " << cumsum << std::endl;
					surface = 0;
					//file.close();

					level++;
					//file.open(infile + std::to_string(level) + "_2.obj");

					//file << str;
					//file << std::endl;
				}
				// now export the node
				//exportBox(*node2.first, file, id++);
				// add the surface 
				h = (*node2.first).bounding_box_proxy().to_bounding_box().surface();
				surface += h;
				cumsum += h;

				// enqueue children, if they exist (!!)
				// if not leaf, add children to the queue
				auto first_child = node2.first->first_child_or_primitive;
				if (node2.first->is_leaf) {
					//typename Bvh::IndexType i;
					//for (i = 0; i < node2.first->primitive_count; i++) {
					//	exportBox(bvh.cnodes[first_child + i], file, id++);
					//}
				}
				else {
					auto* left_child = &bvh.cnodes[first_child + 0];
					auto* right_child = &bvh.cnodes[first_child + 1];

					if (left_child != NULL)
						q2.push(std::make_pair(left_child, node2.second + 1));
					if (right_child != NULL)
						q2.push(std::make_pair(right_child, node2.second + 1));
				}
			}
			statfile.close();
			std::ofstream bigstat;
			bigstat.open(infile + ".txt", std::ios_base::out | std::ios_base::app);
			bigstat << level << " " << cumsum << std::endl; // write the level count and the cumulative surface
			bigstat.close();
		}

		/// Traverse hierarchy and export each level separately
		/// expects a cylinder hierarchy with the root being a box
		void traverseExport(/*std::string filename = infile*/) {
			std::queue<std::pair<Element, int>> q;
			//const auto* node = bvh.cnodes.get();
			int level = 0;
			q.push(std::make_pair(bvh.cnodes.get(), level));
			auto node = q.front();

			// open files
			//std::ofstream bigstat;
			std::ofstream file;
			file.open(infile + "_" + std::to_string(level) + ".obj");
			// copy original .obj
			//std::string str((std::istreambuf_iterator<char>(original)), std::istreambuf_iterator<char>());
			//file << str;
			//file << std::endl;

			int id = 0;
			// count the overall surface area on each level
			Scalar surface = 0;
			Scalar cumsum = 0;
			// and write into a separate file
			std::ofstream statfile;
			statfile.open("stats_cyl.txt");
			Scalar h;
			while (!q.empty()) {
				// pop the first
				node = q.front();
				q.pop();
				// do something to the element
				// while the element level is equal to 'level', write to open file,
				// when it changes, close it and open a new file
				if (node.second != level) {
					// when changing the level, write the summed surface area into the statfile
					statfile << level << " " << (surface) << " " << cumsum << std::endl;
					surface = 0;

					//file.close();
					level++;
					//file.open(infile + "_" + std::to_string((level)) + ".obj");

					//file << str;
					//file << std::endl;
				}
				// now export the node
				//exportBox(*node.first, file, id++);
				// add the surface
				h = (*node.first).bounding_box_proxy().to_bounding_box().surface();
				surface += h;
				cumsum += h;

				// enqueue children, if they exist
				auto first_child = node.first->first_child_or_primitive;
				if (node.first->is_leaf) {
					// if leaf, export every ptimitive's bounding volume
					//typename Bvh::IndexType i;
					//for (i = 0; i < node.first->primitive_count; i++) {
					//	exportBox(bvh.cnodes[first_child + i], file, id++);
					//}
				}
				else {
					auto* left_child = &bvh.cnodes[first_child + 0];
					auto* right_child = &bvh.cnodes[first_child + 1];

					if (left_child != NULL)
						q.push(std::make_pair(left_child, node.second + 1));
					if (right_child != NULL)
						q.push(std::make_pair(right_child, node.second + 1));
				}
			}
			statfile.close();

			std::ofstream bigstat;
			bigstat.open(infile + ".txt", std::ios_base::out | std::ios_base::app);
			bigstat << level << " " << cumsum << std::endl; // write the level count and the cumulative surface
			bigstat.close();
		}

		void traverseExportBox(/*std::string filename = infile*/) {
			std::queue<std::pair<Elem, int>> q;
			int level = 0;
			q.push(std::make_pair(bvh.nodes.get(), level));
			auto node = q.front();

			// open files
			//std::ifstream original;
			std::ofstream file;
			//original.open(infile + ".obj");
			//file.open(infile + "_" + std::to_string(level) + ".obj");
			// copy original .obj
			//std::string str((std::istreambuf_iterator<char>(original)), std::istreambuf_iterator<char>());
			//file << str;
			//file << std::endl;

			int id = 0;

			// count the overall surface area on each level
			Scalar surface = 0;
			Scalar cumsum = 0;
			// and write into a separate file
			std::ofstream statfile;
			statfile.open("stats_box.txt");

			Scalar h;
			while (!q.empty()) {
				// pop the first
				node = q.front();
				q.pop();
				// do something to the element
				// while the element level is equal to level, write to open file,
				// when it changes, close it and open a new file
				if (node.second != level) {
					// when changing the level, write the summed surface area into the statfile
					statfile << level << " " << (surface) << " " << cumsum << std::endl;
					surface = 0;

					//file.close();
					level++;
					//file.open(infile + std::to_string((level)) + ".obj");

					//file << str;
					//file << std::endl;
				}
				// now export the node
				//exportBox(*node.first, file, id++);
				// add the surface 
				h = (*node.first).bounding_box_proxy().to_bounding_box().surface();
				surface += h;
				cumsum += h;

				// enqueue children, if they exist (!!)
				// if not leaf, add children to the queue
				auto first_child = node.first->first_child_or_primitive;
				if (node.first->is_leaf) {
					//typename Bvh::IndexType i;
					//for (i = 0; i < node.first->primitive_count; i++) {
					//	exportBox(bvh.nodes[first_child + i], file, id++);
					//}
				}
				else {
					auto* left_child = &bvh.nodes[first_child + 0];
					auto* right_child = &bvh.nodes[first_child + 1];

					if (left_child != NULL)
						q.push(std::make_pair(left_child, node.second + 1));
					if (right_child != NULL)
						q.push(std::make_pair(right_child, node.second + 1));
				}
			}
			statfile.close();
			std::ofstream bigstat;
			bigstat.open(infile + ".txt", std::ios_base::out | std::ios_base::app);
			bigstat << level << " " << cumsum << std::endl; // write the level count and the cumulative surface
			bigstat.close();
		}
	};
} // namespace bvh

#endif