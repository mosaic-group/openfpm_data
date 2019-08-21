//
// Created by tommaso on 21/8/19.
//

#ifndef OPENFPM_PDATA_PERFORMANCEPLOTS_HPP
#define OPENFPM_PDATA_PERFORMANCEPLOTS_HPP

#include <c++/5/set>
#include <util/common.hpp>

extern char * test_dir;

struct report_sparse_grid_tests
{
    boost::property_tree::ptree graphs;
};

void plotDense2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter);

void plotDense2DZ(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter);

void plotDense3D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
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

void plotInsertSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                              unsigned int &plotCounter);

void plotInsertBlock2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                        unsigned int &plotCounter);

bool isTestInSet(std::set<std::string> &testSet, std::string name)
{
    auto it = testSet.find(name);
    bool res = it != testSet.end();
    std::cout << "isTestInSet - key=\""+ name +"\", found=" << res << std::endl;
    return res;
}

void write_test_report(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet)
{
    const char *perfResultsXmlFile = "SparseGridGpu_performance.xml";

    unsigned int plotCounter = 0;

    plotDense2D(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DZ(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense3D(report_sparsegrid_funcs, testSet, plotCounter);
    plotSparse2DComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotSparse3D(report_sparsegrid_funcs, testSet, plotCounter);
    plotDenseSparse2DComparison(report_sparsegrid_funcs, testSet, plotCounter);

    plotDense2DStencilInsert(report_sparsegrid_funcs, testSet, plotCounter);


    //	report_sparsegrid_funcs.graphs.put("graphs.graph(1).type","line");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).title","SparseGridGPU insert blocked performance");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).x.title","size");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.title","Milion inserts");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.data(0).source","performance.SparseGridGpu(#).insert.Minsert.mean");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).x.data(0).source","performance.SparseGridGpu(#).insert.gridSize.x");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(1).y.data(0).title","line");
//
//	report_sparsegrid_funcs.graphs.put("graphs.graph(2).type","line");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).title","SparseGridGPU insert single performance");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).x.title","size");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.title","Milion inserts");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.data(0).source","performance.SparseGridGpu(#).insertSingle.Minsert.mean");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).x.data(0).source","performance.SparseGridGpu(#).insertSingle.gridSize.x");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(2).y.data(0).title","line");
//
//	report_sparsegrid_funcs.graphs.put("graphs.graph(3).type","line");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).title","SparseGridGPU insert single performance");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).x.title","size");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.title","Milion inserts");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.data(0).source","performance.SparseGridGpu(#).insertStencil.GElems.mean");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).x.data(0).source","performance.SparseGridGpu(#).insertStencil.gridSize.x");
//	report_sparsegrid_funcs.graphs.add("graphs.graph(3).y.data(0).title","line");

    std::string file_xml_results(test_dir);
    file_xml_results += std::string("/") + std::string(perfResultsXmlFile);

    boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
    boost::property_tree::write_xml(file_xml_results, report_sparsegrid_funcs.graphs, std::locale(), settings);

    std::string file_xml_ref(test_dir);
//	file_xml_ref += std::string("/openfpm_pdata/SparseGridGpu_performance_ref.xml"); // This the normal setup
    file_xml_ref += std::string("/") + std::string(perfResultsXmlFile); // This is our setup to get the stdDev on plots

    GoogleChart cg;

    StandardXMLPerformanceGraph(file_xml_results, file_xml_ref, cg, 1);

    addUpdtateTime(cg,1);
    cg.write("SparseGridGpu_performance.html");
}

void
plotDenseSparse2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                            unsigned int plotCounter)
{// 2D Dense sparse comparisons

    std::string dim = "2";
    std::string pattern = "sparse";
    std::string linMode = "N";
    std::string baseDense = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";
    std::string baseSparse = "performance.SparseGridGpu.device.stencil." + pattern + "." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, baseDense + ".gridScaling")
        && isTestInSet(testSet, baseSparse + ".05.gridScaling")
        && isTestInSet(testSet, baseSparse + ".08.gridScaling")
        && isTestInSet(testSet, baseSparse + ".09.gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                baseSparse + ".05.gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, blockEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           baseDense + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           baseSparse + ".05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           baseSparse + ".08.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           baseSparse + ".09.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseDense + ".gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseSparse + ".05.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           baseSparse + ".08.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           baseSparse + ".09.gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "Dense");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "0.5 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "0.8 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title",
                                           "0.9 sparse");
        ++plotCounter;
    }
    if( isTestInSet(testSet, baseDense + ".blockScaling")
        && isTestInSet(testSet, baseSparse + ".05.blockScaling")
        && isTestInSet(testSet, baseSparse + ".08.blockScaling")
        && isTestInSet(testSet, baseSparse + ".09.blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                baseSparse + ".05.blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, gridEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           baseDense + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           baseSparse + ".05.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           baseSparse + ".08.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           baseSparse + ".09.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseDense + ".blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseSparse + ".05.blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           baseSparse + ".08.blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           baseSparse + ".09.blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "Dense");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "0.5 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "0.8 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title",
                                           "0.9 sparse");
        ++plotCounter;
    }
}

void plotSparse3D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter)
{// Sparse 3D
    std::string dim = "3";
    std::string pattern = "sparse";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencil." + pattern + "." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".05.gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + "05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + "05.gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + "05.gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".05.blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil "+pattern+" "+linMode+" "+dim+"D"
                                           +" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + "05.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + "05.blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + "05.blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotSparse2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                            unsigned int &plotCounter)
{// Sparse 2D
    std::string dim = "2";
    std::string pattern = "sparse";
    std::string linMode = "N";
    std::string baseSparse = "performance.SparseGridGpu.device.stencil." + pattern + "." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, baseSparse + ".05.gridScaling")
        && isTestInSet(testSet, baseSparse + ".08.gridScaling")
        && isTestInSet(testSet, baseSparse + ".09.gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                baseSparse + ".05.gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, blockEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           baseSparse + ".05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           baseSparse + ".08.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           baseSparse + ".09.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseSparse + ".05.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseSparse + ".08.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           baseSparse + ".09.gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "0.5 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "0.8 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "0.9 sparse");
        ++plotCounter;
    }
    if( isTestInSet(testSet, baseSparse + ".05.blockScaling")
        && isTestInSet(testSet, baseSparse + ".08.blockScaling")
        && isTestInSet(testSet, baseSparse + ".09.blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                baseSparse + ".05.blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, gridEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           baseSparse + ".05.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           baseSparse + ".08.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           baseSparse + ".09.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseSparse + ".05.blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseSparse + ".08.blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           baseSparse + ".09.blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "0.5 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "0.8 sparse");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "0.9 sparse");
        ++plotCounter;
    }
}

void plotDense3D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 3D
    std::string dim = "3";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil "+linMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2DZ(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                  unsigned int &plotCounter)
{// Dense 2D Z-morton
    std::string dim = "2";
    std::string linMode = "Z";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil "+linMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil "+linMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2DStencilInsert(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencilInsert.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil-insert " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil-insert "+linMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "GFlops");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotInsertSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                        unsigned int &plotCounter)
{
    std::string dim = "2";
    std::string insertMode = "single";
    std::string base = "performance.SparseGridGpu.device.insert.dense." + insertMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert " + insertMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "MInsert/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert "+insertMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "MInsert/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotInsertBlock2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                       unsigned int &plotCounter)
{
    std::string dim = "2";
    std::string insertMode = "block";
    std::string base = "performance.SparseGridGpu.device.insert.dense." + insertMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert " + insertMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "MInsert/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert "+insertMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "MInsert/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
//    report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

#endif //OPENFPM_PDATA_PERFORMANCEPLOTS_HPP
