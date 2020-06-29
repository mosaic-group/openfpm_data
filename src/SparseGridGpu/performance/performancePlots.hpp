//
// Created by tommaso on 21/8/19.
//

#ifndef OPENFPM_PDATA_PERFORMANCEPLOTS_HPP
#define OPENFPM_PDATA_PERFORMANCEPLOTS_HPP

#include <set>
#include <util/common.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <iostream>
#include "Plot/GoogleChart.hpp"
#include "util/performance/performance_util.hpp"

extern char * test_dir;

struct report_sparse_grid_tests
{
    boost::property_tree::ptree graphs;
};

bool isTestInSet(std::set<std::string> &testSet, std::string name);

void write_test_report(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet);

void plotDense2DHost(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotSparse2DHost(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DGetComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DSkeletonComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DZ(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter);

void plotDense2DZComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter);

void plotDense3D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense3DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotSparse2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                            unsigned int &plotCounter);

void plotSparse3D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter);

void
plotDenseSparse2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                            unsigned int plotCounter);

void plotDense2DStencilInsert(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DStencilInsertComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DStencilInplaceInsertComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DStencilInplaceInsertComparison16(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotInsertSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                              unsigned int &plotCounter);

void plotInsertBlock2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                        unsigned int &plotCounter);

void plotGetSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                              unsigned int &plotCounter);

void plotGetNeighbourhood2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                              unsigned int &plotCounter);

#endif //OPENFPM_PDATA_PERFORMANCEPLOTS_HPP
