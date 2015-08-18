/*
 * AdaptiveCellList_test.hpp
 *
 *  Created on: May 4, 2015
 *      Author: i-bird
 */

#ifndef ADAPTIVECELLLIST_TEST_HPP_
#define ADAPTIVECELLLIST_TEST_HPP_

#include "AdaptiveCellList.hpp"
#include "Grid/grid_sm.hpp"

#include <typeinfo>
#include <unordered_set>


// http://stackoverflow.com/questions/20590656/error-for-hash-function-of-pair-of-ints
struct pairhash {
public:
  template <typename T, typename U>
  std::size_t operator()(const std::pair<T, U> &x) const
  {
    return (3 * std::hash<T>()(x.first)) ^ std::hash<U>()(x.second);
  }
};

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
std::unordered_multiset<std::pair<size_t, size_t>, pairhash> allInteractions;

template<int dim, class CellS>
void insertIntoCl(CellS &cl, std::vector<Point<dim, float>> &points) {
	assert(dim >= 2);
	int i = 0;
	for(auto point : points) {
		cl.add(point, i);
		++i;
	}
}

AdaptiveCellList<2, float> getSampleARCL() {
	SpaceBox<2,float> box({0.0f,0.0f},{1.0f,1.0f});
	Point<2,float> org({0.0,0.0});
	AdaptiveCellList<2, float> arcl(box,org);
	
	insertIntoCl<3, AdaptiveCellList<2, float>>(arcl, samplepoints3);
	arcl.construct();
	
	return arcl;
}

BOOST_AUTO_TEST_SUITE( ARCL_and_ARCLvsCL )

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
	timestamp_t c0 = get_timestamp();
	
	auto maxradiusiterator = std::max_element(samplepoints3.begin(), samplepoints3.end(), [](Point<3, float>& a, Point<3, float>& b){return a[2] < b[2];});
	size_t gridsize = std::floor(1.0f / (*maxradiusiterator)[2]);
	//std::cout << "Please choose a gridsize < " << 1.0f / (*maxradiusiterator)[2] << ": " << gridsize << std::endl;
	
	const float epsilon = 0.0001f; // just a little bit bigger, so 1.0 is still inside.
	//const float epsilon = 0.0f; // Iterating over padding cells crashes... due to missing boundary checks, I guess?
	SpaceBox<2,float> box({0.0f,0.0f},{1.0f+epsilon,1.0f+epsilon});
	size_t div[2] = {gridsize,gridsize};
	grid_sm<2,void> g_info(div);
	Point<2,float> org({0.0,0.0});
	const float pad = 1;
	CellList<2, float> cl(box,div,org,pad, 30000);
	
	insertIntoCl<2, CellList<2, float>>(cl, samplepoints2); // shares indices with ~3, so just use the 2d one for a simple celllist
	
	timestamp_t c1 = get_timestamp();
	
	std::cout << "Creation time (us): " << c1-c0 << std::endl;
	
	
	size_t interactions_count = 0, interaction_candidates = 0;
	
	size_t cell, i1, i2;
	
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
						//std::cout << i1 << " and " << i2 << " with p2 = " << samplepoints3[i2].toString() << std::endl;
						allInteractions.insert(std::make_pair(i1, i2));
						++interactions_count;
					}
					++iter_inner;
				}
				//std::cout << std::endl;
				++iter;
			}
		}
	}
	
	timestamp_t t1 = get_timestamp();
	std::cout << "Interactions time (us): " << t1-t0 << std::endl;
	std::cout << "Found interactions: " << interactions_count
			<< " of " << interaction_candidates
			<< " candidates (" << (static_cast<float>(interaction_candidates) / interactions_count)
			<< "x)." << std::endl;
}

BOOST_AUTO_TEST_CASE( find_inserted_items_in_arcl)
{
	auto arcl = getSampleARCL();
	
	for(auto& p: samplepoints3) {
		size_t cellindex = arcl.findCellIndex(p);
		
		// check here, whether the calculated cell really contains our point!
		bool found = false;
		
		auto iter = arcl.getCellContents(cellindex);
		for(auto childiter = iter.first; childiter != iter.second; ++childiter)
			if(childiter->first == p)
				found = true;
		
		BOOST_REQUIRE(found);
	}
	
	//auto printresult = arcl.printTree(0,0,0);
	//std::cout << printresult.first << printresult.second << std::endl;
}

BOOST_AUTO_TEST_CASE( findcellindex_and_findcellcenter_combine)
{
	auto arcl = getSampleARCL();
	
	size_t max = 164193736414550; // = theoretically we should choose -1, which is 18446744073709551615 on my system, but on high levels float fails, so we stop earlier.
	size_t sqrt_max = sqrt(max);
	//max = (1l << 32) + 1l;
	for(size_t i = 0; i < 10000000; i++)
		BOOST_REQUIRE_EQUAL(arcl.findCellIndex(arcl.findCellCenter(i)), i);
	std::cout << "Checked until " << 10000000 << std::endl;
	for(size_t i = 100000000; i < sqrt_max+10000; i += 777777)
		BOOST_REQUIRE_EQUAL(arcl.findCellIndex(arcl.findCellCenter(i)), i);
	std::cout << "Checked until " << sqrt_max+100000 << std::endl;
	for(size_t i = sqrt_max+10000; i < max; i += sqrt_max)
		BOOST_CHECK_EQUAL(arcl.findCellIndex(arcl.findCellCenter(i)), i);
	std::cout << "Checked until " << max << std::endl;
}

BOOST_AUTO_TEST_CASE( get_all_interactions_arcl)
{
	timestamp_t c0 = get_timestamp();
	auto arcl = getSampleARCL();
	timestamp_t c1 = get_timestamp();
	std::cout << "Creation time (us): " << c1-c0 << std::endl;
	
	size_t interactions_count = 0, interaction_candidates = 0;
	
	size_t i1, i2;
	
	timestamp_t t0 = get_timestamp();
	
	///*
	for(std::pair<Point<3,float>, size_t>& p1 : arcl) {
		i1 = p1.second;
		auto iter_inner = arcl.getNNIterator<FAST>(p1.first); //full interactions
		while (iter_inner.isNext()) {
			auto p2 = iter_inner.get();
			i2 = p2.second;
			++interaction_candidates;
			if(i1 != i2 && is_in_radius(p1.first, p2.first)) {
				++interactions_count;
				auto it = allInteractions.find(std::make_pair(i1, i2));
				BOOST_REQUIRE(it != allInteractions.end());
				allInteractions.erase(it);
				//std::cout << i1 << " and " << i2 << " with p2 = " << p2.first.toString() << std::endl;
			} //else std::cout << i1 << " | | " << i2 << " with p2 = " << p2.first.toString() << std::endl;
			++iter_inner;
		}
		//std::cout << std::endl;
	}
	BOOST_REQUIRE_EQUAL(allInteractions.size(), 0);
	//*/
	
	//std::cout << "0.252 ^ 2:" << std::endl;
	//arcl.getNNIterator<FAST>(Point<3,float>({0.02f,0.02f,0.05f}));
	
	//for(size_t i = 0; i < 10; i++)
	//	arcl.getNNIterator<FAST>(arcl.findCellCenter(i));
	
	
	timestamp_t t1 = get_timestamp();
	std::cout << "Interaction time (us): " << t1-t0 << std::endl;
	std::cout << "Found interactions: " << interactions_count
			<< " of " << interaction_candidates
			<< " candidates (" << (static_cast<float>(interaction_candidates) / interactions_count)
			<< "x)." << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* ADAPTIVECELLLIST_TEST_HPP_ */
