/*
 * GoogleChart.hpp
 *
 *  Created on: Jan 9, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PLOT_GOOGLECHART_HPP_
#define OPENFPM_DATA_SRC_PLOT_GOOGLECHART_HPP_

#include <fstream>

#define GGRAPH_COLUMS 1

/*! \brief Google chart options
 *
 */
struct GCoptions
{
	//! Title of the chart
	std::string title;
	//! Y axis name
	std::string yAxis;
	//! X axis name
	std::string xAxis;

	//! Type of chart (list of the option can be founded in Google Chart API for seriesType)
	//! Possible options are:
	//!  'line', 'area', 'bars', 'candlesticks', and 'steppedArea'
	//! default: line
	std::string stype;

	//! Extended series options
	//! Example {5: {type: 'line'}} specify that the series number 5 must be represented
	//! with a line
	std::string stypeext;

	size_t width=900;
	size_t heigh=500;

	GCoptions & operator=(const GCoptions & opt)
	{
		return *this;
	}
};

struct GGraph
{
	// TypeOfGraph
	size_t type;

	// data
	std::string data;

	// option
	std::string option;

	// Google chart option
	GCoptions opt;
};

/////////////////// Constants strings usefull to construct the HTML page //////////

const std::string begin_data ="<html>\n\
  <head>\n\
    <script type=\"text/javascript\" src=\"https://www.gstatic.com/charts/loader.js\"></script>\n\
    <script type=\"text/javascript\">\n\
      google.charts.load('current', {'packages':['corechart']});\n\
      google.charts.setOnLoadCallback(drawVisualization);\n\
\n\
\n\
      function drawVisualization() {\n";

const std::string end_data="]);\n\n";

const std::string begin_div = "}</script>\n\
</head>\n\
<body>\n";

const std::string div_end = "</body>\n\
</html>\n";

/////////////////////////////////////////////////////////////////////

/*! It convert an array y or a set of two array x,y into a Google chart
 *
 */
class GoogleChart
{
	// set of graphs
	openfpm::vector<GGraph> set_of_graphs;

	// set inject HTML;
	openfpm::vector<std::string> injectHTML;

	/*! \brief Given X and Y vector
	 *
	 *
	 */
	template<typename X, typename Y, typename Yn> std::string get_data(const openfpm::vector<X> & x, const openfpm::vector<Y> & y, const openfpm::vector<Yn> & yn, const GCoptions & opt)
	{
		std::stringstream data;

		// we require that the number of x elements are the same as y elements

		if (x.size() != y.size())
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " vector x and the vector y must have the same number of elements " << x.size() << "!=" << y.size() << "\n";

		data << "['";
		data << opt.xAxis << "'";
		for (size_t i = 0 ; i < yn.size() ; i++)
			data << ",'" << std::to_string(yn.get(i)) << "'";

		data << "],\n";

		// Write the values
		for (size_t i = 0 ; i < y.size() ; i++)
		{
			data << "[";
			data << "'" << std::to_string(x.get(i)) << "'";
			for (size_t j = 0 ; j < y.get(i).size() ; j++)
				data << "," << std::to_string(y.get(i).get(j));

			data << "],\n";
		}

		return data.str();
	}

	std::string get_option(const GCoptions & opt)
	{
		std::stringstream str;
		str << "title : '" << opt.title << "',\n";
	    str << "vAxis: {title: '" << opt.yAxis <<  "'},\n";
	    str << "hAxis: {title: '" << opt.xAxis << "'},\n";
	    str << "seriesType: '" << opt.stype << "',\n";
	    if (opt.stypeext.size() != 0)
	    	str << "series: " << opt.stypeext << "\n";

	    return str.str();
	}

	/*! \brief Add a graph data variable
	 *
	 * \param of file out
	 * \param i id
	 * \param data string
	 *
	 */
	void addData(std::ofstream & of, size_t i, const std::string & data)
	{
		of << "var data";
		of << i;
		of << " = google.visualization.arrayToDataTable([\n";
		of << data;
		of << "]);\n";
	}

	/*! \brief Add an option data variable
	 *
	 * \param of file out
	 * \param i id
	 * \param opt string
	 *
	 */
	void addOption(std::ofstream & of, size_t i, const std::string & opt)
	{
		of << "var options";
		of << i;
		of << "= {\n";
		of << opt;
		of << "};\n";
	}

	/*! \brief Add a draw div section
	 *
	 * \param of file out
	 * \param i id
	 *
	 */
	void addDrawDiv(std::ofstream & of, size_t i)
	{
		of << "var chart = new google.visualization.ComboChart(document.getElementById('chart_div";
		of << i;
		of << "'));chart.draw(data";
		of << i;
		of << ", options";
		of << i;
		of << ");\n";
	}

	/*! \brief Add a div section
	 *
	 * \param i id
	 * \param gc GoogleChart option
	 *
	 */
	void addDiv(std::ofstream & of, size_t i, const GCoptions & gc)
	{
		of << "<div id=\"chart_div";
		of << i;
		of << "\" style=\"width: ";
		of << gc.width;
		of << "px; height: ";
		of << gc.heigh;
		of << "px;\"></div>\n";
	}

public:

	GoogleChart()
	{
		injectHTML.add();
	}

	/*! \brief Add a colums graph
	 *
	 * \param y A vector of vector of values (numbers) the size of y indicate how many columns
	 *          has the graph while the internal vector can store multiple datasets
	 *
	 */
	template<typename Y> void AddColumsGraph(openfpm::vector<Y> & y)
	{
		openfpm::vector<std::string> x;
		x.resize(y.size());

		AddColumsGraph<std::string,Y>(x,y);
	}

	/*! \brief Add a colums graph
	 *
	 * \param y A vector of vector of values (numbers) the size of y indicate how many columns
	 *          has the graph while the internal vector can store multiple datasets
	 *
	 * \param x Give a name or number to each colums, so can be a string or a number
	 *
	 */
	template<typename X, typename Y> void AddColumsGraph(openfpm::vector<X> & x, openfpm::vector<Y> & y)
	{
		GCoptions opt;

		openfpm::vector<std::string> yn;

		if (y.size() != 0)
			yn.resize(y.get(0).size());

		AddColumsGraph<X,Y,std::string>(x,y,yn,opt);
	}

	/*! \brief Add a colums graph
	 *
	 * \param y A vector of vector of values (numbers) the size of y indicate how many columns
	 *          has the graph while the internal vector can store multiple datasets
	 *
	 * \param x Give a name or number to each colums, so can be a string or a number
	 *
	 * \param yn Give a name to each dataset
	 *
	 */
	template<typename X, typename Y, typename Yn> void AddColumsGraph(openfpm::vector<X> & x, openfpm::vector<Y> & y, openfpm::vector<Yn> & yn)
	{
		GCoptions opt;

		AddColumsGraph(x,y,yn,opt);
	}

	/*! \brief Add a colums graph
	 *
	 * \param y A vector of vector of values (numbers) the size of y indicate how many columns
	 *          has the graph while the internal vector can store multiple datasets
	 *
	 * \param x Give a name or number to each colums, so can be a string or a number
	 *
	 * \param yn Give a name to each dataset
	 *
	 * \param opt Graph options
	 *
	 */
	template<typename X, typename Y, typename Yn> void AddColumsGraph(openfpm::vector<X> & x, openfpm::vector<Y> & y, openfpm::vector<Yn> & yn , const GCoptions & opt)
	{
		set_of_graphs.add();
		injectHTML.add();

		// Check that all the internal vector has the same number of elements

		if (y.size() != 0)
		{
			size_t sz = y.get(0).size();
			for (size_t i = 0; i < y.size() ; i++)
			{
				if (y.get(i).size() != sz)
					std::cerr << __FILE__ << ":" << __LINE__ << " error all the elements in the y vector must have the same numbers, element " << i << ": " << y.get(i).size() << " " << " mismatch the numbers of elements at 0: " << sz << "/n";
			}
		}

		set_of_graphs.last().type = GGRAPH_COLUMS;
		set_of_graphs.last().data = get_data(x,y,yn,opt);
		set_of_graphs.last().option = get_option(opt);
		set_of_graphs.last().opt = opt;
	}

	/*! \brief Add HTML text
	 *
	 * \param html add html text in the page
	 *
	 */
	void addHTML(const std::string & html)
	{
		injectHTML.last() = html;
	}

	/*! \brief It write the html file
	 *
	 * \param file output html file
	 *
	 */
	void write(std::string file)
	{
		// Open a file

		std::ofstream of(file);

		// Check if the file is open
		if (of.is_open() == false)
		{std::cerr << "Error cannot create the HTML file: " + file + "\n";}

		// write the file

		of << begin_data;

		for (size_t i = 0 ; i < set_of_graphs.size() ; i++)
			addData(of,i,set_of_graphs.get(i).data);

		for (size_t i = 0 ; i < set_of_graphs.size() ; i++)
			addOption(of,i,set_of_graphs.get(i).option);

		for (size_t i = 0 ; i < set_of_graphs.size() ; i++)
			addDrawDiv(of,i);

		of << begin_div;

		of << injectHTML.get(0);

		for (size_t i = 0 ; i < set_of_graphs.size() ; i++)
		{
			addDiv(of,i,set_of_graphs.get(i).opt);
			of << injectHTML.get(i+1);
		}

		of << div_end;

		of.close();
	}
};

#endif /* OPENFPM_DATA_SRC_PLOT_GOOGLECHART_HPP_ */
