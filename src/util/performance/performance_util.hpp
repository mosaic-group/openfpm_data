/*
 * performance_util.hpp
 *
 *  Created on: Jul 21, 2019
 *      Author: i-bird
 */

#ifndef PERFORMANCE_UTIL_HPP_
#define PERFORMANCE_UTIL_HPP_

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

static void addUpdtateTime(GoogleChart & cg, int np)
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

    std::stringstream str;

    std::string commit;
    commit = exec("git rev-parse HEAD");

    str << "<h3>Updated: " << now->tm_mday << "/" << now->tm_mon + 1 << "/" << now->tm_year+1900 << "     " << now->tm_hour << ":" << now->tm_min << ":"
    		               << now->tm_sec << "  commit: " << commit << "   run with: " << np << " processes" << std::endl;

	cg.addHTML(str.str());
}

static void createCommitFile(std::string & tmp)
{
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );

    std::stringstream str;

    std::string commit;
    commit = exec("git rev-parse HEAD");

	ofstream f("comit_f_" + tmp);
	f << commit << std::endl;

	f.close();
}

static inline void warning_set(int & warning_level, double mean, double mean_ref, double sigma)
{
	int warning_level_candidate;

	if (mean - mean_ref < -2.0*sigma )
		warning_level_candidate = -1;
	else if (mean - mean_ref < 2.0*sigma)
		warning_level_candidate = 0;
	else if (mean - mean_ref < 3.0*sigma)
		warning_level_candidate = 1;
	else
		warning_level_candidate = 2;

	if (warning_level_candidate > warning_level)
		warning_level = warning_level_candidate;
}

static inline void addchartarea(std::string & chart_area, int lvl)
{
	std::string color;

	if (lvl == -1)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#00FF00',\
		        strokeWidth: 6\
		    }\
		}");
	}
	else if (lvl == 0)
	{
		// NOTHING TO DO
	}
	else if (lvl == 1)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#FFFF00',\
		        strokeWidth: 6\
		    }\
		}");
	}
	else if (lvl == 2)
	{
		chart_area = std::string(",chartArea: {\
		    backgroundColor: {\
		        stroke: '#FF0000',\
		        strokeWidth: 6\
		    }\
		}");
	}

}

/*! \brief Check of a string is a float
 *
 * \param something strign to check if it is a number
 *
 */
static bool isFloat(const std::string &someString)
{
	using boost::lexical_cast;
	using boost::bad_lexical_cast;

	try
	{boost::lexical_cast<float>(someString);}
	catch (bad_lexical_cast &)
	{return false;}

	return true;
}

static void StandardXMLPerformanceGraph(std::string file_xml,
		                      std::string file_xml_ref,
							  GoogleChart & cg,
							  const double deviationMultiplier = 3.0)
{
    // Create empty property tree object
	boost::property_tree::ptree tree_measure;

    // Parse the XML into the property tree.
	boost::property_tree::read_xml(file_xml, tree_measure);

    // Create empty property tree object
    boost::property_tree::ptree tree_reference;

    // We check if exist the reference file. If does not exist copy the file_xml into the report folder
    if (boost::filesystem::exists(file_xml_ref) == false )
    {
    	boost::filesystem::copy_file(file_xml,file_xml_ref);
    }

    // Parse the JSON into the property tree.
    boost::property_tree::read_xml(file_xml_ref, tree_reference);

    // First we check for graphs

//    try
    {
    	boost::property_tree::ptree childs = tree_measure.get_child("graphs");

    	for (auto & c: childs)
    	{
    		std::string type = c.second.template get<std::string>("type","");
    		// discorver the number of points
    		int number = 0;
    		openfpm::vector<std::string> yn;

    		while (1)
    		{
    			if (c.second.template get<std::string>("y.data(" + std::to_string(number) + ").title","") == "")
    			{break;}
    			yn.add(c.second.template get<std::string>("y.data(" + std::to_string(number) + ").title",
    			        "line" + std::to_string(number)));
    			yn.add("interval");
    			yn.add("interval");
    			number++;
    		}

    		bool is_log_x = c.second.template get<bool>("options.log_x",false);
    		bool is_log_y = c.second.template get<bool>("options.log_y",false);

    		// We process the graph
    		std::string title = c.second.template get<std::string>("title","");
    		std::string x_title = c.second.template get<std::string>("x.title","");
    		std::string y_title = c.second.template get<std::string>("y.title","");

    		// This is optional

    	    std::string preHTML = c.second.template get<std::string>("preHTML","");
    		std::string postHTML = c.second.template get<std::string>("postHTML","");

    		openfpm::vector<openfpm::vector<double>> x;
    		openfpm::vector<openfpm::vector<std::string>> xs;
    		openfpm::vector<openfpm::vector<double>> y;

    		openfpm::vector<openfpm::vector<double>> x_ref;
    		openfpm::vector<openfpm::vector<double>> y_ref_up;
    		openfpm::vector<openfpm::vector<double>> y_ref_dw;

    		int warning_level = -1;
    		bool is_literal_x = false;

    		for (size_t i = 0 ; i < number ; i++)
    		{
    			x.add();
    			xs.add();
    			y.add();

    			x_ref.add();
    			y_ref_up.add();
    			y_ref_dw.add();

    			std::string xv = c.second.template get<std::string>("x.data(" + std::to_string(i) + ").source","");
    			std::string yv = c.second.template get<std::string>("y.data(" + std::to_string(i) + ").source","");

    			// Get the numbers

    			int j = 0;
    			while (1)
    			{
    				std::string xv_ = xv;
    				std::string yv_ = yv;

    				int xpos = xv_.find("#");
    				int ypos = yv_.find("#");

    				xv_.replace(xpos,1,std::to_string(j));
    				yv_.replace(ypos,1,std::to_string(j));

    				if (tree_measure.template get<std::string>(xv_,"") == "")
    				{
    					if (j == 0)
    					{
    						std::cout << "WARNING: Not found " << xv_ << " in file: " << file_xml << std::endl;
    					}
    					break;
    				}

    				std::string tmp = tree_measure.template get<std::string>(xv_,"");

    				double x_val;
    				std::string x_val_str;
    				if (isFloat(tmp) == true)
    				{x_val = tree_measure.template get<double>(xv_,0.0);}
    				else
    				{
    					is_literal_x = true;
    					xs.last().add(tmp);
    				}

    				double y_val = tree_measure.template get<double>(yv_,0.0);

    				double x_val_ref = tree_reference.template get<double>(xv_,0.0);
    				double y_val_ref = tree_reference.template get<double>(yv_,0.0);

    				if (y_val_ref == 0.0)
    				{
    					std::cout << "WARNING: " << yv_ << " does not exist on the reference file: " << file_xml_ref << std::endl;
    				}

    				ypos = yv_.find(".mean");
    				std::string yv_dev = yv_.replace(ypos,5,".dev");
    				double y_val_dev_ref = tree_reference.template get<double>(yv_,0.0);

    				if (y_val_dev_ref == 0.0)
    				{
    					std::cout << "WARNING: " << yv_ << " does not exist on the reference file: " << file_xml_ref << std::endl;
    				}

    				x.last().add(x_val);
    				y.last().add(y_val);

    				x_ref.last().add(x_val_ref);
                    y_ref_dw.last().add(y_val_ref - deviationMultiplier * y_val_dev_ref);
    				y_ref_up.last().add(y_val_ref + deviationMultiplier * y_val_dev_ref);

					warning_set(warning_level,y_val,y_val_ref,y_val_dev_ref);

    				j++;
    			}
    		}

    		if (type == "line")
    		{
    			GCoptions opt;

    			std::string chart_area;
    			addchartarea(chart_area,warning_level);
    			opt.curveType = c.second.template get<std::string>("interpolation","function");
//    			opt.curveType = c.second.template get<std::string>("interpolation","none");

    			if (is_log_x == true)
    			{
    				opt.more = GC_X_LOG + "," + GC_ZOOM + chart_area;
    			}
    			else
        		{
        			opt.more = GC_ZOOM + chart_area;
        		}

                if (is_log_y == true)
                {
                    opt.more = GC_Y_LOG + "," + GC_ZOOM + chart_area;
                }
                else
                {
                    opt.more = GC_ZOOM + chart_area;
                }

    			opt.title = title;
    			opt.xAxis = x_title;
    			opt.yAxis = y_title;

    			if (is_literal_x == false)
    			{
					if (x.size() == 1)
					{
							cg.AddLines(yn,opt,x.get(0),y.get(0),
											x_ref.get(0),y_ref_dw.get(0),
											x_ref.get(0),y_ref_up.get(0));
					}
					if (x.size() == 2)
					{
							cg.AddLines(yn,opt,x.get(0),y.get(0),
										   x_ref.get(0),y_ref_dw.get(0),x_ref.get(0),y_ref_up.get(0),
										   x.get(1),y.get(1),
										   x_ref.get(1),y_ref_dw.get(1),x_ref.get(1),y_ref_up.get(1));
					}
					if (x.size() == 3)
					{
							cg.AddLines(yn,opt,x.get(0),y.get(0),
										   x_ref.get(0),y_ref_dw.get(0),x_ref.get(0),y_ref_up.get(0),
										   x.get(1),y.get(1),
										   x_ref.get(1),y_ref_dw.get(1),x_ref.get(1),y_ref_up.get(1),
										   x.get(2),y.get(2),
										   x_ref.get(2),y_ref_dw.get(2),x_ref.get(2),y_ref_up.get(2));
					}
					if (x.size() == 4)
					{
							cg.AddLines(yn,opt,x.get(0),y.get(0),
											   x_ref.get(0),y_ref_dw.get(0),x_ref.get(0),y_ref_up.get(0),
											   x.get(1),y.get(1),
											   x_ref.get(1),y_ref_dw.get(1),x_ref.get(1),y_ref_up.get(1),
											   x.get(2),y.get(2),
											   x_ref.get(2),y_ref_dw.get(2),x_ref.get(2),y_ref_up.get(2),
											   x.get(3),y.get(3),
											   x_ref.get(3),y_ref_dw.get(3),x_ref.get(3),y_ref_up.get(3));
					}
					if (x.size() == 5)
					{
						cg.AddLines(yn,opt,x.get(0),y.get(0),
										   x_ref.get(0),y_ref_dw.get(0),x_ref.get(0),y_ref_up.get(0),
										   x.get(1),y.get(1),
										   x_ref.get(1),y_ref_dw.get(1),x_ref.get(1),y_ref_up.get(1),
										   x.get(2),y.get(2),
										   x_ref.get(2),y_ref_dw.get(2),x_ref.get(2),y_ref_up.get(2),
										   x.get(3),y.get(3),
										   x_ref.get(3),y_ref_dw.get(3),x_ref.get(3),y_ref_up.get(3),
										   x.get(4),y.get(4),
										   x_ref.get(4),y_ref_dw.get(4),x_ref.get(4),y_ref_up.get(4));
					}
					if (x.size() == 6)
					{
						cg.AddLines(yn,opt,x.get(0),y.get(0),
										   x_ref.get(0),y_ref_dw.get(0),x_ref.get(0),y_ref_up.get(0),
										   x.get(1),y.get(1),
										   x_ref.get(1),y_ref_dw.get(1),x_ref.get(1),y_ref_up.get(1),
										   x.get(2),y.get(2),
										   x_ref.get(2),y_ref_dw.get(2),x_ref.get(2),y_ref_up.get(2),
										   x.get(3),y.get(3),
										   x_ref.get(3),y_ref_dw.get(3),x_ref.get(3),y_ref_up.get(3),
										   x.get(4),y.get(4),
										   x_ref.get(4),y_ref_dw.get(4),x_ref.get(4),y_ref_up.get(4),
										   x.get(5),y.get(5),
										   x_ref.get(5),y_ref_dw.get(5),x_ref.get(5),y_ref_up.get(5));
					}
    			}
				else
				{
					openfpm::vector<openfpm::vector<double>> y_tmp;

					y_tmp.resize(x.get(0).size());

					for (size_t i = 0 ; i < x.size() ; i++)
					{
						for (size_t j = 0 ; j < x.get(i).size() ; j++)
						{
							y_tmp.get(j).add(y.get(i).get(j));
							y_tmp.get(j).add(y_ref_dw.get(i).get(j));
							y_tmp.get(j).add(y_ref_up.get(i).get(j));
						}
					}

					cg.AddLinesGraph(xs.get(0),y_tmp,yn,opt);
				}
    		}
    	}
    }
//    catch (std::exception e)
//    {
//    	std::cout << __FILE__ << ":" << __LINE__ << " Error: invalid xml for performance test " << e.what() << std::endl;
//    }

}


#endif /* PERFORMANCE_UTIL_HPP_ */
