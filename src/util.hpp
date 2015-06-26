/*
 * util.hpp
 *
 *  Created on: May 7, 2015
 *      Author: Pietro Incardona
 */

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <boost/iostreams/device/mapped_file.hpp>


/*! \brief Compare two files, return true if they match
 *
 * \param file1 path1
 * \param file2 path2
 *
 * \return true if they match
 *
 */
bool compare(std::string file1, std::string file2)
{
    boost::iostreams::mapped_file_source f1(file1);
    boost::iostreams::mapped_file_source f2(file2);

    if( f1.size() == f2.size() && std::equal(f1.data(), f1.data() + f1.size(), f2.data()) )
        return true;
    else
        return false;
}

struct RGB
{
	float R;
	float G;
	float B;

	std::string toString()
	{
		return std::to_string(R) + " " + std::to_string(G) + " " + std::to_string(B);
	}
};

/*! \brief Return the color sampled from one group
 *
 * groups:
 *
 * 0: Red
 * 1: Green
 * 2: Blue
 * 3: Yellow
 * 4: Cyan
 * 5: Magenta
 * 6: Orange
 * 7: Chartreuse-Green
 * 8: Spring Green
 * 9: Azure
 * 10: Violet
 * 11: Rose
 *
 * \param group
 *
 */

struct RGB getColor(int group, std::uniform_real_distribution<float> & d, std::default_random_engine & g)
{
	struct RGB col;

	if (group == 0)
	{
		float s = d(g);
		col.R = s/2 + 0.5;
		col.G = 0.0;
		col.B = 0.0;
	}
	else if (group == 1)
	{
		float s = d(g);
		col.R = 0.0;
		col.G = s/2 + 0.5;
		col.B = 0.0;
	}
	else if (group == 2)
	{
		float s = d(g);
		col.R = 0.0;
		col.G = 0.0;
		col.B = s;
	}
	else if (group == 3)
	{
		float s = d(g);
		col.R = s/2 + 0.5;
		col.G = s/2 + 0.5;
		col.B = 0.0;
	}
	else if (group == 4)
	{
		float s = d(g);
		col.R = s/2 + 0.5;
		col.G = 0.0;
		col.B = s/2 + 0.5;
	}
	else if (group == 5)
	{
		float s = d(g);
		col.R = 0.0;
		col.G = s/2 + 0.5;
		col.B = s/2 + 0.5;
	}
	else if (group == 6)
	{
		float s = d(g);
		col.R = s/2 + 0.5;
		col.G = s/4 + 0.5;
		col.B = 0.0;
	}
	else if (group == 7)
	{
		float s = d(g);
		col.R = s/4 + 0.5;
		col.G = s/2 + 0.5;
		col.B = 0.0;
	}
	else if (group == 8)
	{
		float s = d(g);
		col.R = 0.0;
		col.G = s/2 + 0.5;
		col.B = s/4 + 0.5;
	}
	else if (group == 9)
	{
		float s = d(g);
		col.R = 0.0;
		col.G = s/4 + 0.5;
		col.B = s/2 + 0.5;
	}
	else if (group == 10)
	{
		float s = d(g);
		col.R = s/4 + 0.5;
		col.G = 0.0;
		col.B = s/2 + 0.5;
	}
	else if (group == 11)
	{
		float s = d(g);
		col.R = s/2 + 0.5;
		col.G = 0.0;
		col.B = s/4 + 0.5;
	}

	return col;
}

/*! \brief Check if one string end with a particular string
 *
 * \param fullString string to check
 * \param ending ending string to check
 *
 */
bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length())
    {return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));}
    else
    {return false;}
}

#endif /* UTIL_HPP_ */
