/*
 * performancePlots.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: i-bird
 */

#include "performancePlots.hpp"

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

    plotDense2DHost(report_sparsegrid_funcs, testSet, plotCounter);
    plotSparse2DHost(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2D(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DZ(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DZComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DGetComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DSkeletonComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense3D(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense3DComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotSparse2DComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotSparse3D(report_sparsegrid_funcs, testSet, plotCounter);
    plotDenseSparse2DComparison(report_sparsegrid_funcs, testSet, plotCounter);

    plotDense2DStencilInsert(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DStencilInsertComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DStencilInplaceInsertComparison(report_sparsegrid_funcs, testSet, plotCounter);
    plotDense2DStencilInplaceInsertComparison16(report_sparsegrid_funcs, testSet, plotCounter);
    plotInsertSingle2D(report_sparsegrid_funcs, testSet, plotCounter);
    plotInsertBlock2D(report_sparsegrid_funcs, testSet, plotCounter);

    plotGetSingle2D(report_sparsegrid_funcs, testSet, plotCounter);
    plotGetNeighbourhood2D(report_sparsegrid_funcs, testSet, plotCounter);



    std::string file_xml_results(test_dir);
    file_xml_results += std::string("/") + std::string(perfResultsXmlFile);

    boost::property_tree::xml_writer_settings<std::string> settings(' ', 4);
    boost::property_tree::write_xml(file_xml_results, report_sparsegrid_funcs.graphs, std::locale(), settings);

    std::string file_xml_ref(test_dir);
//	file_xml_ref += std::string("/openfpm_pdata/SparseGridGpu_performance_ref.xml"); // This the normal setup
    file_xml_ref += std::string("/") + std::string(perfResultsXmlFile); // This is our setup to get the stdDev on plots

    GoogleChart cg;

    StandardXMLPerformanceGraph(file_xml_results, file_xml_ref, cg, 1);

    addUpdateTime(cg,1,"data","SparseGridGpu_performance");
    cg.write("SparseGridGpu_performance.html");
}

void plotDenseSparse2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                            unsigned int plotCounter)
{// 2D Dense sparse comparisons

    std::string dim = "2";
    std::string pattern = "sparse";
    std::string linMode = "N";
    std::string baseDense = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";
    std::string baseSparse = "performance.SparseGridGpu.device.stencil." + pattern + "." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, baseDense + ".8.gridScaling")
        && isTestInSet(testSet, baseSparse + ".05.gridScaling")
        && isTestInSet(testSet, baseSparse + ".08.gridScaling")
        && isTestInSet(testSet, baseSparse + ".09.gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                baseSparse + ".05.gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil dense-" + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, blockEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           baseDense + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           baseSparse + ".05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           baseSparse + ".08.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           baseSparse + ".09.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseDense + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseSparse + ".05.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           baseSparse + ".08.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           baseSparse + ".09.gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

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
                                           "SparseGridGPU stencil dense-" + pattern + " " + linMode + " " + dim + "D"
                                           + " grid scaling performance, gridEdge=" + std::to_string(bes));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
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
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + "05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + "05.gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + "05.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + "05.blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
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
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
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
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2DHost(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.host.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}
void plotSparse2DHost(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.host.stencil.sparse." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".05.gridScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".05.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".05.gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".05.gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=" + std::to_string(bes));
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".05.blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil "+linMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".05.blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".05.blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".05.blockScaling(0).gridSize.x"));
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".4.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                   base + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                   base + ".16.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                   base + ".32.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");


        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title", "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title", "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title", "blockEdge=32");

        ++plotCounter;
    }
}

void plotDense2DSkeletonComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.skeleton.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU skeleton stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Elem/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".4.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                   base + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                   base + ".16.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                   base + ".32.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");


        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title", "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title", "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title", "blockEdge=32");

        ++plotCounter;
    }
}

void plotDense2DGetComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencilGet.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU GET stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".4.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                   base + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                   base + ".16.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                   base + ".32.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");


        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title", "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title", "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title", "blockEdge=32");

        ++plotCounter;
    }
}

void plotDense2DZComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "Z";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".4.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                   base + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                   base + ".16.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                   base + ".32.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");


        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title", "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title", "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title", "blockEdge=32");

        ++plotCounter;
    }
}

void plotDense3DComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "3";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".4.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                   base + ".8.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");


        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title", "blockEdge=8");

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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotDense2DStencilInsertComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string base = "performance.SparseGridGpu.device.stencilInsert.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil-insert " + linMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                   base + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   base + ".16.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".16.gridScaling(#).gridSize.x");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "blockEdge=16");

        ++plotCounter;
    }
}

void plotDense2DStencilInplaceInsertComparison(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string baseInplace = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";
    std::string baseInsert = "performance.SparseGridGpu.device.stencilInsert.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, baseInplace + ".8.gridScaling")
        && isTestInSet(testSet, baseInsert + ".8.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil in-place vs. insert " + linMode + " " + dim + "D" +
                                           " - blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                   baseInplace + ".8.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   baseInsert + ".8.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseInplace + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseInsert + ".8.gridScaling(#).gridSize.x");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "In-place");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "Insert");

        ++plotCounter;
    }
}

void plotDense2DStencilInplaceInsertComparison16(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                 unsigned int &plotCounter)
{// Dense 2D
    std::string dim = "2";
    std::string linMode = "N";
    std::string baseInplace = "performance.SparseGridGpu.device.stencil.dense." + linMode + "." + dim + "D";
    std::string baseInsert = "performance.SparseGridGpu.device.stencilInsert.dense." + linMode + "." + dim + "D";

    if( isTestInSet(testSet, baseInplace + ".16.gridScaling")
        && isTestInSet(testSet, baseInsert + ".16.gridScaling")
      )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU stencil in-place vs. insert " + linMode + " " + dim + "D" +
                                           " - blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G flops/s");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                   baseInplace + ".16.gridScaling(#).GFlops.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                   baseInsert + ".16.gridScaling(#).GFlops.mean");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           baseInplace + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           baseInsert + ".16.gridScaling(#).gridSize.x");

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);

        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title", "In-place");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title", "Insert");

        ++plotCounter;
    }
}

void plotInsertSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                        unsigned int &plotCounter)
{
    std::string dim = "2";
    std::string insertMode = "single";
    std::string base = "performance.SparseGridGpu.device.insert.dense." + insertMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert " + insertMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "M Inserts/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           base + ".4.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           base + ".8.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
//        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
//                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "blockEdge=8");
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "M Inserts/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
   report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
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

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
            )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU insert " + insertMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "M Inserts/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           base + ".4.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           base + ".8.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           base + ".16.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                           base + ".32.gridScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
//        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
//                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title",
                                           "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title",
                                           "blockEdge=32");
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
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "M Inserts/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Minsert.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotGetSingle2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                       unsigned int &plotCounter)
{
    std::string dim = "2";
    std::string getMode = "single";
    std::string base = "performance.SparseGridGpu.device.get.dense." + getMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
            )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU get " + getMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Gets/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           base + ".4.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           base + ".8.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           base + ".16.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                           base + ".32.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
//        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
//                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title",
                                           "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title",
                                           "blockEdge=32");
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU get "+getMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Gets/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

void plotGetNeighbourhood2D(report_sparse_grid_tests &report_sparsegrid_funcs, std::set<std::string> &testSet,
                       unsigned int &plotCounter)
{
    std::string dim = "2";
    std::string getMode = "neighbourhood";
    std::string base = "performance.SparseGridGpu.device.get.dense." + getMode + "." + dim + "D";

    if( isTestInSet(testSet, base + ".2.gridScaling")
        && isTestInSet(testSet, base + ".4.gridScaling")
        && isTestInSet(testSet, base + ".8.gridScaling")
        && isTestInSet(testSet, base + ".16.gridScaling")
        && isTestInSet(testSet, base + ".32.gridScaling")
            )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU get " + getMode + " " + dim + "D" +
                                           " grid scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title", "GridEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Gets/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".2.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).source",
                                           base + ".4.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).source",
                                           base + ".8.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).source",
                                           base + ".16.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).source",
                                           base + ".32.gridScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".2.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(1).source",
                                           base + ".4.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(2).source",
                                           base + ".8.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(3).source",
                                           base + ".16.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(4).source",
                                           base + ".32.gridScaling(#).gridSize.x");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x", true);
//        int bes = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
//                base + ".gridScaling(0).blockSize"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "blockEdge=2");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(1).title",
                                           "blockEdge=4");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(2).title",
                                           "blockEdge=8");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(3).title",
                                           "blockEdge=16");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(4).title",
                                           "blockEdge=32");
        ++plotCounter;
    }
    if( isTestInSet(testSet, base + ".blockScaling") )
    {
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").type", "line");
        report_sparsegrid_funcs.graphs.put("graphs.graph(" + std::to_string(plotCounter) + ").interpolation", "none");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").title",
                                           "SparseGridGPU get "+getMode+" "+dim+"D"+" block scaling performance");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.title",
                                           "BlockEdgeSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.title", "G Gets/s");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).source",
                                           base + ".blockScaling(#).Gget.mean");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").x.data(0).source",
                                           base + ".blockScaling(#).blockSize");
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").options.log_x",true);
        int ges = static_cast<int>( report_sparsegrid_funcs.graphs.template get<double>(
                base + ".blockScaling(0).gridSize.x"));
        report_sparsegrid_funcs.graphs.add("graphs.graph(" + std::to_string(plotCounter) + ").y.data(0).title",
                                           "gridEdge=" + std::to_string(ges));
        ++plotCounter;
    }
}

