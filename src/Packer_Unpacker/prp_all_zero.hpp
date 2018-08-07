/*
 * prp_all_zero.hpp
 *
 *  Created on: Mar 3, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PACKER_UNPACKER_PRP_ALL_ZERO_HPP_
#define OPENFPM_DATA_SRC_PACKER_UNPACKER_PRP_ALL_ZERO_HPP_

#include "util/variadic_to_vmpl.hpp"

//! Structure to convert a variadic template into boost::mpl::vector
template<typename T, bool is_zero, int ... prp>
struct prp_all_zero
{
	//! mpl vector
	typedef typename to_boost_vmpl<prp...>::type type;
};

//! Structure to convert a variadic template into boost::mpl::vector
//! when sizeof...(prp) == 0
template<typename T, int ... prp>
struct prp_all_zero<T,true,prp ... >
{
	//! mpl vector
	typedef typename boost::mpl::range_c<int,0, boost::mpl::size<typename T::type>::type::value >::type type;
};


#endif /* OPENFPM_DATA_SRC_PACKER_UNPACKER_PRP_ALL_ZERO_HPP_ */
