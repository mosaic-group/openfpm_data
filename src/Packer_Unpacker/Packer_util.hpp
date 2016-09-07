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
	template<typename T, typename Memory=HeapMemory, typename layout=typename memory_traits_lin<T>::type, template<typename> class layout_base=memory_traits_lin , typename grow_p=grow_policy_double, unsigned int impl=vect_isel<T>::value> class vector;
}

///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// FUNCTORS FOR ENCAP ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////


//A functor for call_aggregatePackRequest
template<typename encap, typename Mem>
struct call_packRequest_enc_functor
{

	encap & obj;
	size_t & req;

	call_packRequest_enc_functor(encap & obj, size_t & req)
	:obj(obj),req(req)
	{}

	//! It calls the pack request for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename encap::T_type,T>::type obj_type;

		Packer<obj_type,Mem>::packRequest(obj.template get<T::value>(),req);
	}
};


//Calls a pack request in nested way
template<typename encap, typename Mem, int ... prp>
struct call_encapPackRequest
{
	static inline void call_packRequest(encap & obj, size_t & req)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<encap,sizeof...(prp) == 0,prp...>::type b_prp;

		call_packRequest_enc_functor<encap,Mem> functor(obj,req);

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

template<typename obj_type, typename Mem>
struct call_packRequest_agg_functor
{

	const obj_type & obj;
	size_t & req;

	call_packRequest_agg_functor(const obj_type & obj, size_t & req)
	:obj(obj),req(req)
	{}

	//! It calls the pack request for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename obj_type::type,T>::type obj_t;

		//for (size_t i = 0; i < obj_type::max_prop ; i++)
		Packer<obj_t,Mem>::packRequest(obj.template get<T::value>(),req);
	}
};


//Calls a pack request in nested way
template<typename obj_type, typename Mem, int ... prp>
struct call_aggregatePackRequest
{
	static inline void call_packRequest(const obj_type & obj, size_t & req)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<obj_type,sizeof...(prp) == 0,prp...>::type b_prp;

		call_packRequest_agg_functor<obj_type,Mem> functor(obj,req);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};

//A functor for call_aggregatePack
template<typename obj_type, typename Mem>
struct call_pack_agg_functor
{
	ExtPreAlloc<Mem> & mem;
	const obj_type & obj;
	Pack_stat & sts;

	call_pack_agg_functor(ExtPreAlloc<Mem> & mem, const obj_type & obj, Pack_stat & sts)
	:mem(mem), obj(obj), sts(sts)
	{
	}

	//! It calls the packer for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename obj_type::type,T>::type obj_t;

		//for (size_t i = 0; i < obj.size(); i++)
			Packer<obj_t,Mem>::pack(mem,obj.template get<T::value>(),sts);
	}
};

//Calls a packer in nested way
template<typename obj_type, typename Mem, int ... prp>
struct call_aggregatePack
{
	static inline void call_pack(const obj_type & obj, ExtPreAlloc<Mem> & mem, Pack_stat & sts)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<obj_type,sizeof...(prp) == 0,prp...>::type b_prp;

		call_pack_agg_functor<obj_type,Mem> functor(mem,obj,sts);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};


//A functor for call_aggregateUnpack
template<typename obj_type, typename Mem>
struct call_unpack_agg_functor
{
	ExtPreAlloc<Mem> & mem;
	const obj_type & obj;
	Unpack_stat & ps;

	call_unpack_agg_functor(ExtPreAlloc<Mem> & mem, const obj_type & obj, Unpack_stat & ps)
	:mem(mem), obj(obj), ps(ps)
	{
	}

	//! It calls the unpacker for each property
	template<typename T>
	inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<typename obj_type::type,T>::type obj_t;

		//for (size_t i = 0; i < obj.size(); i++)
			Unpacker<obj_t,Mem>::unpack(mem,obj.template get<T::value>(),ps);
	}
};

//Calls an unpacker in nested way
template<typename obj_type, typename Mem, int ... prp>
struct call_aggregateUnpack
{
	static inline void call_unpack(const obj_type & obj, ExtPreAlloc<Mem> & mem, Unpack_stat & ps)
	{
		//Property sequence into boost::mpl::range_c or boost::mpl::vector, depending on sizeof...(prp)
		typedef typename prp_all_zero<obj_type,sizeof...(prp) == 0,prp...>::type b_prp;

		call_unpack_agg_functor<obj_type,Mem> functor(mem,obj,ps);

		//Apply functor for each property
		boost::mpl::for_each_ref<b_prp>(functor);
	}
};

#endif /* PACKER_UTIL_HPP_ */
