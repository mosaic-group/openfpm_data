#define BOOST_DISABLE_ASSERTS

#include "config.h"
#define BOOST_TEST_MODULE "C++ test module for OpenFPM_data project"
#include <boost/test/included/unit_test.hpp>

#include <iostream>
#include <boost/mpl/int.hpp>
#include <typeinfo>

// Include tests

/*
#include "timer_util_test.hpp"
#include "Grid/grid_key_dx_expression_unit_tests.hpp"
#include "Point_test_unit_tests.hpp"
#include "util/util_test.hpp"
#include "Space/SpaceBox_unit_tests.hpp"
#include "Space/Shape/Box_unit_tests.hpp"
#include "NN/CellList/CellList_test.hpp"
#include "Vector/vector_unit_tests.hpp"
#include "hypercube_unit_test.hpp"
#include "Graph/graph_unit_tests.hpp"
#include "Grid/grid_unit_tests.hpp"
*/
#include "NN/AdaptiveCellList/AdaptiveCellList_test.hpp"


BOOST_AUTO_TEST_SUITE( CellListComparison )

#include <sys/time.h>
typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp ()
{
	struct timeval now;
	gettimeofday (&now, NULL);
	return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

#include <fstream>
template<int dim>
std::vector<Point<dim, float>> getSamplePoints() {
	std::vector<Point<dim, float>> result;
	std::vector<short> img(1080*1080*3, 0); // 1080x1080 rgb with black bg
	result.reserve(61650);
	
	std::ifstream rsource, xsource, ysource;
	rsource.open("../../sampledata/r.txt", std::ifstream::in);
	xsource.open("../../sampledata/x.txt", std::ifstream::in);
	ysource.open("../../sampledata/y.txt", std::ifstream::in);
	
	std::string rline, xline, yline;
	
	for(; std::getline(rsource, rline); )
	{
		std::getline(xsource, xline);
		std::getline(ysource, yline);
		
		std::istringstream rin(rline);
		std::istringstream xin(xline);
		std::istringstream yin(yline);

		float r, x, y;
		
		rin >> r;
		xin >> x;
		yin >> y;
		
		if (dim == 2)
			result.push_back({x, y});
		else
			result.push_back({y,x, r}); // TODO: change this back ;)
		
		float size = 0.05f;
		int cr,cg,cb, index, nx,ny;
		cr = 255.0f*(r/0.09f);
		cg = 255.0f*(1.0f-r/0.09f);
		cb = 100;
		int bx = static_cast<int>(x*1079.0f);
		int by = static_cast<int>(y*1079.0f);
		for(int dy = -r*1079.0f; dy < r*1079.0f; dy++)
		for(int dx = -r*1079.0f; dx < r*1079.0f; dx++) {
			nx = bx + dx;
			ny = by + dy;
			if(nx < 0 || nx > 1079 || ny < 0 || ny > 1079) continue;
			if(dx*dx+dy*dy > size*r*1079.0f*r*1079.0f) continue;
			index = ny * 1080 + nx;
			img[3*index] = cr;
			img[3*index+1] = cg;
			img[3*index+2] = cb-cb*(dx*dx+dy*dy)/(r*1079.0f*r*1079.0f);
		}
	}
	
	rsource.close();
	xsource.close();
	ysource.close();
	
	std::ofstream imgfile;
	imgfile.open("../../sampledata/img.ppm", std::ofstream::out);
	imgfile << "P3\n1080 1080\n255\n";
	int index;
	for(int y=0;y<1080;y++) {
		for(int x=0;x<1080;x++) {
			index = y*1080+x;
			imgfile << img[3*index] << " " << img[3*index+1] << " " << img[3*index+2] << " ";
		}
		imgfile << "\n";
	}
	imgfile.close();
	
	return result;
	
	//return std::vector<Point<dim, float>>({
	//	Point<dim, float>({0.3f, 0.06f, 0.0656f})
	//});
}


std::vector<Point<2, float>> samplepoints2 = getSamplePoints<2>();
std::vector<Point<3, float>> samplepoints3 = getSamplePoints<3>();

template<int dim, class CellS>
void insertIntoCl(CellS &cl, std::vector<Point<dim, float>> &points) {
	assert(dim >= 2);
	int i = 0;
	for(auto point : points) {
		cl.add(point, i);
		++i;
	}
}

BOOST_AUTO_TEST_CASE( celllist_realloc)
{
	std::cout << "2D Cell List:" << std::endl;
	
	//Space where is living the Cell list
	SpaceBox<2,float> box({0.0f,0.0f},{1.0f,1.0f});

	// Subdivisions
	size_t div[2] = {16,16};
	
	// grid info
	grid_sm<2,void> g_info(div);

	// Origin
	Point<2,float> org({0.0,0.0});

	// Cell lists (slot numbers are based on div = 16^2)
	CellList<2, float> cl1(box,div,org,1); // fallback slot = 16
	CellList<2, float> cl2(box,div,org,1, 6716); // one realloc needed
	CellList<2, float> cl3(box,div,org,1, 13431); // precise number
	CellList<2, float> cl4(box,div,org,1, 100000); // wild overestimate


	CellList<2, float> cl0(box,div,org,1);
	for(float i=0.0f; i < 1.0; i += 0.05f)for(float j=0.0f; j < 1.0; j += 0.05f)
	cl0.add(Point<2, float>({i,j}), 4711);
	//cl0.add(Point<2, float>({1.0,1.0}), 4711);
	//cl0.add(Point<2, float>({1.3,0.1}), 4711);
	//cl0.add(Point<2, float>({2,2}), 4711);

	// insert points
	timestamp_t t0 = get_timestamp();
	insertIntoCl<2, CellList<2, float>>(cl1, samplepoints2);
	timestamp_t t1 = get_timestamp();
	insertIntoCl<2, CellList<2, float>>(cl2, samplepoints2);
	timestamp_t t2 = get_timestamp();
	insertIntoCl<2, CellList<2, float>>(cl3, samplepoints2);
	timestamp_t t3 = get_timestamp();
	insertIntoCl<2, CellList<2, float>>(cl4, samplepoints2);
	timestamp_t t4 = get_timestamp();
	
	int sum1 = 0;
	for (unsigned int i=0; i<div[0]+2; i++) for (unsigned int j=0; j<div[1]+2; j++) sum1 += cl1.getNelements(cl1.getGrid().LinId(grid_key_dx<2>(i, j)));
	int sum2 = 0;
	for (unsigned int i=0; i<div[0]+2; i++) for (unsigned int j=0; j<div[1]+2; j++) sum2 += cl2.getNelements(cl2.getGrid().LinId(grid_key_dx<2>(i, j)));
	int sum3 = 0;
	for (unsigned int i=0; i<div[0]+2; i++) for (unsigned int j=0; j<div[1]+2; j++) sum3 += cl3.getNelements(cl3.getGrid().LinId(grid_key_dx<2>(i, j)));
	int sum4 = 0;
	for (unsigned int i=0; i<div[0]+2; i++) for (unsigned int j=0; j<div[1]+2; j++) sum4 += cl4.getNelements(cl4.getGrid().LinId(grid_key_dx<2>(i, j)));

	/*
	for (int i=0; i<div[0]+2; i++) {
		for (int j=0; j<div[1]+2; j++) {
			::printf("%6d  ",cl0.getNelements(cl0.getGrid().LinId(grid_key_dx<2>(i, j))));
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	
	for (int i=0; i<div[0]+2; i++) {
		for (int j=0; j<div[1]+2; j++) {
			::printf("%6d  ",cl1.getNelements(cl1.getGrid().LinId(grid_key_dx<2>(i, j))));
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	size_t s = samplepoints2.size();
	
	std::cout << "Inserted " << s << " elements" << std::endl;
	
	BOOST_REQUIRE_EQUAL(s, sum1);
	BOOST_REQUIRE_EQUAL(s, sum2);
	BOOST_REQUIRE_EQUAL(s, sum3);
	BOOST_REQUIRE_EQUAL(s, sum4);
	
	std::cout << "10 reallocs (us): " << t1-t0 << std::endl;
	std::cout << "1 realloc (us): " << t2-t1 << std::endl;
	std::cout << "Without realloc, precise (us): " << t3-t2 << std::endl;
	std::cout << "Without realloc, overestimate (us): " << t4-t3 << std::endl;
}

inline bool is_in_radius(const Point<3, float> &p1, const Point<3, float> &p2) {
	Point<3, float> diff = p1 - p2;
	//std::cout << p1.toString() << "  -  " << p2.toString() << "  =  " << diff.toString() << std::endl << std::flush;
	//std::cout << p2.get(2) << std::endl << std::flush;
	//std::cout << diff[0]*diff[0] + diff[1]*diff[1] <<std::flush<< " <= " << (p1[2] < p2[2] ? p1[2] : p2[2]) << std::endl << std::flush;
	float minradius = (p1[2] < p2[2] ? p1[2] : p2[2]);
	return diff[0]*diff[0] + diff[1]*diff[1] <= minradius*minradius;
}

BOOST_AUTO_TEST_CASE( get_all_interactions_classiccelllist)
{
	return;
	
	auto maxradiusiterator = std::max_element(samplepoints3.begin(), samplepoints3.end(), [](Point<3, float>& a, Point<3, float>& b){return a[2] < b[2];});
	size_t gridsize = std::floor(1.0f / (*maxradiusiterator)[2]);
	//gridsize = 7;
	std::cout << "Please choose a gridsize < " << 1.0f / (*maxradiusiterator)[2] << ": " << gridsize << std::endl;
	
	const float epsilon = 0.0001f; // just a little bit bigger, so 1.0 is still inside.
	//const float epsilon = 0.0f; // Iterating over padding cells crashes... due to missing boundary checks, I guess?
	SpaceBox<2,float> box({0.0f,0.0f},{1.0f+epsilon,1.0f+epsilon});
	size_t div[2] = {gridsize,gridsize};
	grid_sm<2,void> g_info(div);
	Point<2,float> org({0.0,0.0});
	const float pad = 1;
	CellList<2, float> cl(box,div,org,pad, 30000);
	
	insertIntoCl<2, CellList<2, float>>(cl, samplepoints2); // shares indices with ~3, so just use the 2d one for a simple celllist
	
	size_t interactions_count = 0, interaction_candidates = 0;
	
	size_t cell, i1, i2;
	
	std::cout << "ready lets go" << std::endl;
	timestamp_t t0 = get_timestamp();
	
	for (int i=pad; i<div[0]+pad; i++) {
		for (int j=pad; j<div[1]+pad; j++) {
			cell = cl.getGrid().LinId(grid_key_dx<2>(i, j));
			//std::cout << "getting iterator for cell " << cell << " with elems: " << nels << std::endl;
			auto iter = cl.getNNIterator<FAST>(cell); //full interactions
			while (iter.isNext()) {
				i2 = iter.get();
				//std::cout << "Neighbor found: " << i2 << std::endl;
				CellIterator<CellList<2, float>> iter_inner(cell, cl);
				while (iter_inner.isNext()) {
					i1 = iter_inner.get();
					++interaction_candidates;
					if(i1 != i2 && is_in_radius(samplepoints3[i1], samplepoints3[i2])) {
						//std::cout << i1 << " and " << i2 << std::endl;
						++interactions_count;
					}//else std::cout << i1 << " AND " << i2 << std::endl;
					++iter_inner;
				}
				//std::cout << std::endl;
				++iter;
			}
		}
	}
	
	timestamp_t t1 = get_timestamp();
	std::cout << "that was easy (us): " << t1-t0 << std::endl;
	std::cout << "found interactions: " << interactions_count
			<< " of " << interaction_candidates
			<< " candidates (" << (static_cast<float>(interaction_candidates) / interactions_count)
			<< "x)." << std::endl;
}

BOOST_AUTO_TEST_CASE( get_all_interactions_arcelllist)
{
	SpaceBox<2,float> box({0.0f,0.0f},{1.0f,1.0f});
	Point<2,float> org({0.0,0.0});
	AdaptiveCellList<2, float> arcl(box,org);
	
	insertIntoCl<3, AdaptiveCellList<2, float>>(arcl, samplepoints3);
	arcl.construct();
	
	size_t interactions_count = 0, interaction_candidates = 0;
	
	size_t i1, i2;
	
	std::cout << "ready lets go (adaptive)" << std::endl;
	
	timestamp_t t0 = get_timestamp();
	
	auto iter = arcl.getNNIterator<FAST>(4711); //full interactions
	while (iter.isNext()) {
		i2 = iter.get();
		//std::cout << "Neighbor found: " << i2 << std::endl;
		auto iter_inner = arcl.getNNIterator<FAST>(4711); //full interactions
		while (iter_inner.isNext()) {
			i1 = iter_inner.get();
			++interaction_candidates;
			if(i1 != i2 && is_in_radius(samplepoints3[i1], samplepoints3[i2])) {
				++interactions_count;
				//std::cout << i1 << " and " << i2 << std::endl;
			}
			++iter_inner;
		}
		//std::cout << std::endl;
		++iter;
	}
	
	timestamp_t t1 = get_timestamp();
	std::cout << "that was easy (us): " << t1-t0 << std::endl;
	std::cout << "found interactions (should be 475560): " << interactions_count
			<< " of " << interaction_candidates
			<< " candidates (" << (static_cast<float>(interaction_candidates) / interactions_count)
			<< "x)." << std::endl;
	
	auto printresult = arcl.printTree(0,0,0);
	std::cout << printresult.first << printresult.second << std::endl;
	
	for(auto& p: samplepoints3)
		arcl.findCellIndex(p);
	arcl.findCellIndex(Point<3,float>({0.0f,0.0f,0.066f}));
	
	std::cout << "0: ";
	for (auto& i : arcl.findChildrenIndices(0))
		std::cout << i << " ";
	std::cout << std::endl;
	
	std::cout << "1: ";
	for (auto& i : arcl.findChildrenIndices(1))
		std::cout << i << " ";
	std::cout << std::endl;
	
	std::cout << "5: ";
	for (auto& i : arcl.findChildrenIndices(5))
		std::cout << i << " ";
	std::cout << std::endl;
	
	std::cout << "6: ";
	for (auto& i : arcl.findChildrenIndices(6))
		std::cout << i << " ";
	std::cout << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()