/*
 * MetaParser.hpp
 *
 *  Created on: Mar 3, 2019
 *      Author: i-bird
 */

#ifndef METAPARSER_HPP_
#define METAPARSER_HPP_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/bind/bind.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>

typedef boost::program_options::options_description MetaParser_options;
namespace MetaParser_def = boost::program_options;

class MetaParser
{
	boost::program_options::options_description desc;

	boost::program_options::variables_map vm;

public:

	explicit MetaParser(boost::program_options::options_description & desc)
	:desc(desc)
	{
	}

	/*! \brief Parse the string of options
	 *
	 * \param string of options
	 *
	 */
	bool parse(std::string & opts)
	{
		std::istringstream iss(opts);

		// Parse mocked up input.
		boost::program_options::store(boost::program_options::parse_config_file(iss,desc), vm);
		boost::program_options::notify(vm);

		return true;
	}

	/*! \brief Return the option opt in value
	 *
	 * \param opt option to check
	 * \param value where to store the value of the option
	 *
	 * \return true if the option has been set
	 *
	 */
	template<typename T>
	bool getOption(std::string opt,T & value)
	{
		if (vm.count(opt))
		{
			value = vm[opt].as<T>();
			return true;
		}

		return false;
	}

	/*! \brief clear everything you parsed so far
	 *
	 *
	 */
	void clear()
	{
		vm.clear();
	}
};


#endif /* METAPARSER_HPP_ */
