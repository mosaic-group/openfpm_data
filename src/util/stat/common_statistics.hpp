/*
 * common_statistics.hpp
 *
 *  Created on: Jul 20, 2019
 *      Author: i-bird
 */

#ifndef COMMON_STATISTICS_HPP_
#define COMMON_STATISTICS_HPP_


/*! \brief Standard deviation
 *
 * \param measures set of measures
 * \param mean the mean of the measures
 *
 * \return the standard deviation
 *
 */
static inline void standard_deviation(openfpm::vector<double> measures, double & mean, double & dev)
{
	mean = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
		mean += measures.get(i);
	mean /= measures.size();

	dev = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
		dev += (measures.get(i) - mean)*(measures.get(i) - mean);

	dev = sqrt(dev / (measures.size() - 1));
}

/*! \brief Standard deviation
 *
 * \param measures set of measures
 * \param mean the mean of the measures
 *
 * \return the standard deviation
 *
 */
static inline void standard_deviation(std::vector<double> measures, double & mean, double & dev)
{
	mean = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
	{
		mean += measures[i];
	}

	mean /= measures.size();

	dev = 0;
	for (size_t i = 0 ; i < measures.size() ; i++)
		dev += (measures[i] - mean)*(measures[i] - mean);

	dev = sqrt(dev / (measures.size() - 1));
}


#endif /* COMMON_STATISTICS_HPP_ */
