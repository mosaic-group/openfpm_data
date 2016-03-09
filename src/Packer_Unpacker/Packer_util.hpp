/*
 * Packer_util.hpp
 *
 *  Created on: Jan 25, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef PACKER_UTIL_HPP_
#define PACKER_UTIL_HPP_

#include "prp_all_zero.hpp"
#include "Packer_Unpacker/Pack_selector.hpp"
#include "memory/HeapMemory.hpp"
#include "Vector/map_vector_grow_p.hpp"
#include "Vector/vect_isel.hpp"
#include "memory_ly/memory_conf.hpp"

template<typename T, typename Mem, int pack_type=Pack_selector<T>::value > class Packer;
template<typename T, typename Mem, int pack_type=Pack_selector<T>::value > class Unpacker;

namespace openfpm
{
	template<typename T, typename Memory=HeapMemory, typename grow_p=grow_policy_double, unsigned int impl=vect_isel<T>::value> class vector;
}

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// FUNCTORS FOR ENCAP ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


//A functor for call_aggregatePackRequest
template<typename encap, typename Mem>
struct call_packRequest_enc_functor
{

	encap & obj;
	std::vector<size_t> & pap_prp;

	call_packRequest_enc_functor(encap & obj, std::vector<size_t> & pap_prp)
	:obj(obj),pap_prp(pap_prp)
	{}

	//! It calls the pack request for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename encap::T_type,T>::type obj_type;

		Packer<obj_type,Mem>::packRequest(obj.template get<T::value>(),pap_prp);
	}
};


//Calls a pack request in nested way
template<typename encap, typename Mem, int ... prp>
struct call_encapPackRequest
{
	static inline void call_packRequest(encap & obj, std::vector<size_t> & pap_prp)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<encap,sizeof...(prp) == 0,prp...>::type b_prp;

		call_packRequest_enc_functor<encap,Mem> functor(obj,pap_prp);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};

//A functor for call_aggregatePack
template<typename encap, typename Mem>
struct call_pack_enc_functor
{
	ExtPreAlloc<Mem> & mem;
	const encap & obj;
	Pack_stat & sts;

	call_pack_enc_functor(ExtPreAlloc<Mem> & mem,const encap & obj, Pack_stat & sts)
	:mem(mem), obj(obj), sts(sts)
	{
	}

	//! It calls the packer for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename encap::type,T>::type obj_type;

		Packer<obj_type,Mem>::pack(mem,obj.template get<T::value>(),sts);
	}
};

//Calls a packer in nested way
template<typename encap, typename Mem, int ... prp>
struct call_encapPack
{
	static inline void call_pack(const encap & obj, ExtPreAlloc<Mem> & mem, Pack_stat & sts)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<typename encap::T_type,sizeof...(prp) == 0,prp...>::type b_prp;

		call_pack_enc_functor<encap,Mem> functor(mem,obj,sts);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};


//A functor for call_aggregateUnpack
template<typename encap, typename Mem, int ... prp>
struct call_unpack_encap_functor
{
	ExtPreAlloc<Mem> & mem;
	const encap & obj;
	Unpack_stat & ps;

	call_unpack_encap_functor(ExtPreAlloc<Mem> & mem, const encap & obj, Unpack_stat & ps)
	:mem(mem), obj(obj), ps(ps)
	{
	}

	//! It calls the unpacker for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename encap::type,T>::type obj_type;

		Unpacker<obj_type,Mem>::unpack(mem,obj.template get<T::value>(),ps);
	}
};

//Calls an unpacker in nested way
template<typename encap, typename Mem, int ... prp>
struct call_encapUnpack
{
	static inline void call_unpack(encap & obj, ExtPreAlloc<Mem> & mem, Unpack_stat & ps)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<encap,sizeof...(prp) == 0,prp...>::type b_prp;

		call_unpack_encap_functor<encap,Mem,prp ... > functor(mem,obj,ps);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};


////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////// FUNCTORS FOR AGGREGATE ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

//A functor for call_aggregatePackRequest
template<typename Aggr, typename Mem, typename grow_p>
struct call_packRequest_agg_functor
{

	openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj;
	std::vector<size_t> & pap_prp;

	call_packRequest_agg_functor(openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj, std::vector<size_t> & pap_prp)
	:obj(obj),pap_prp(pap_prp)
	{}

	//! It calls the pack request for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename Aggr::type,T>::type obj_type;

		for (size_t i = 0; i < obj.size(); i++)
			Packer<obj_type,Mem>::packRequest(obj.template get<T::value>(i),pap_prp);
	}
};


//Calls a pack request in nested way
template<typename T, typename Mem, typename grow_p, int ... prp>
struct call_aggregatePackRequest
{
	static inline void call_packRequest(openfpm::vector<T,Mem,grow_p,OPENFPM_NATIVE> & obj, std::vector<size_t> & pap_prp)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<T,sizeof...(prp) == 0,prp...>::type b_prp;

		call_packRequest_agg_functor<T,Mem,grow_p> functor(obj,pap_prp);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};

//A functor for call_aggregatePack
template<typename Aggr, typename Mem, typename grow_p>
struct call_pack_agg_functor
{
	ExtPreAlloc<Mem> & mem;
	openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj;
	Pack_stat & sts;

	call_pack_agg_functor(ExtPreAlloc<Mem> & mem, openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj, Pack_stat & sts)
	:mem(mem), obj(obj), sts(sts)
	{
	}

	//! It calls the packer for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename Aggr::type,T>::type obj_type;

		for (size_t i = 0; i < obj.size(); i++)
			Packer<obj_type,Mem>::pack(mem,obj.template get<T::value>(i),sts);
	}
};

//Calls a packer in nested way
template<typename T, typename Mem, typename grow_p, int ... prp>
struct call_aggregatePack
{
	static inline void call_pack(openfpm::vector<T,Mem,grow_p,OPENFPM_NATIVE> & obj, ExtPreAlloc<Mem> & mem, Pack_stat & sts)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<T,sizeof...(prp) == 0,prp...>::type b_prp;

		call_pack_agg_functor<T,Mem,grow_p> functor(mem,obj,sts);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};


//A functor for call_aggregateUnpack
template<typename Aggr, typename Mem, typename grow_p>
struct call_unpack_agg_functor
{
	ExtPreAlloc<Mem> & mem;
	openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj;
	Unpack_stat & ps;

	call_unpack_agg_functor(ExtPreAlloc<Mem> & mem, openfpm::vector<Aggr,Mem,grow_p,OPENFPM_NATIVE> & obj, Unpack_stat & ps)
	:mem(mem), obj(obj), ps(ps)
	{
	}

	//! It calls the unpacker for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename Aggr::type,T>::type obj_type;

		for (size_t i = 0; i < obj.size(); i++)
			Unpacker<obj_type,Mem>::unpack(mem,obj.template get<T::value>(i),ps);
	}
};

//Calls an unpacker in nested way
template<typename T, typename Mem, typename grow_p, int ... prp>
struct call_aggregateUnpack
{
	static inline void call_unpack(openfpm::vector<T,Mem,grow_p,OPENFPM_NATIVE> & obj, ExtPreAlloc<Mem> & mem, Unpack_stat & ps)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<T,sizeof...(prp) == 0,prp...>::type b_prp;

		call_unpack_agg_functor<T,Mem,grow_p> functor(mem,obj,ps);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};

//A functor for has_aggregatePack
template<typename Aggr>
struct has_pack_agg_functor
{
	bool at_least_one = false;

	//! It calls the "has_pack" checking for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename Aggr::type,T>::type obj_type;

		at_least_one |= has_pack<obj_type>::value;
	}
};

//Returns true if at least one property of aggregate has "pack()" member, false is not
template<typename T, int ... prp>
struct has_aggregatePack
{
	static inline bool has_pack()
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<T,sizeof...(prp) == 0,prp...>::type b_prp;

		has_pack_agg_functor<T> functor;

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);

		return functor.at_least_one;
	}
};


#endif /* PACKER_UTIL_HPP_ */
