/*
 * has_pack_encap.hpp
 *
 *  Created on: Mar 3, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_PACK_ENCAP_HPP_
#define OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_PACK_ENCAP_HPP_

#include "has_pack_agg.hpp"
#include "prp_all_zero.hpp"


#include <boost/fusion/mpl.hpp>


template<class T, int ... prp>
struct has_pack_encap
{
	typedef typename prp_all_zero<T,sizeof...(prp) == 0,prp...>::type vprp;
	//typedef typename to_boost_vmpl<prp...>::type vprp;
    typedef typename has_pack_agg_impl<typename T::type,sizeof ... (prp), boost::mpl::bool_<false> , vprp>::result result;
};


#endif /* OPENFPM_DATA_SRC_PACKER_UNPACKER_HAS_PACK_ENCAP_HPP_ */
