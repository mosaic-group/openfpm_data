/*
 * trash_bin.hpp
 *
 *  Created on: Nov 24, 2014
 *      Author: i-bird
 */

#ifndef TRASH_BIN_HPP_
#define TRASH_BIN_HPP_


/*! \brief Edge class that encapsulate an object T
 *
 * Edge class that encapsulate an object T that store the property of the edge
 * It is basically used to append a field size_t
 * to the original object T, this field is the neighborhood vertex in the adjacency list
 *
 */

template<typename T>
class E_p
{
public:
	//! insert a size_t property to the boost::vector
	typedef typename boost::fusion::result_of::push_front<typename T::type, size_t>::type E_a;

	//! type for the edge
	typedef E_a type;

	E_a data;
};

/*! \brief Vertex class that encapsulate an object T
 *
 * Vertex class that encapsulate an object T that store the property of the vertex
 * It is basically used to insert a field size_t
 * to the original object T, this field is the end list point in the adjacency list
 *
 */

template<typename T>
class V_p
{
public:
	//! insert a size_t property to the boost::vector
	typedef typename boost::fusion::result_of::push_front<typename T::type, size_t>::type V_a;

	//! type for the vertex
	typedef V_a type;

	V_a data;
};


#endif /* TRASH_BIN_HPP_ */


/*! \brief Get the reference of the selected element
 *
 * Get the reference of the selected element
 *
 * \param p property to get (is an integer)
 * \param v1 grid_key that identify the element in the grid
 *
 */
/*		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & getBoostVector(grid_key_dx<dim> & v1)
{
#ifdef MEMLEAK_CHECK
	check_valid(&boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1))));
#endif
	return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
}*/


template<unsigned int prp>
struct red_max
{
	template<typename aggr_vect, typename encap_type>
	static void red(aggr_vect & aggr, encap_type & ec)
	{
		auto op1 = ec.template get<prp>();
		auto op2 = aggr.template get<prp>();

		aggr.template get<prp>() = (op2 < op1)?op1:op2;
	}

	template<typename red_type, typename aggr_vect, typename BlockReduce>
	static red_type red_final(red_type * rt, aggr_vect & aggr)
	{
		return BlockReduce(rt).Max(aggr.template get<prp>());
	}
};

template<unsigned int prp>
struct red_min
{
	template<typename aggr_vect, typename encap_type>
	static void red(aggr_vect & aggr, encap_type & ec)
	{
		auto op1 = ec.template get<prp>();
		auto op2 = aggr.template get<prp>();

		aggr.template get<prp>() = (op2 > op1)?op1:op2;
	}


	template<typename red_type, typename aggr_vect, typename BlockReduce>
	static red_type red_final(red_type * rt, aggr_vect & aggr)
	{
		return BlockReduce(rt).Min(aggr.template get<prp>());
	}
};

template<unsigned int prp>
struct red_sum
{
	template<typename aggr_vect, typename encap_type>
	static void red(aggr_vect & aggr, encap_type & ec)
	{
		aggr.template get<prp>() += ec.template get<prp>();
	}


	template<typename red_type, typename aggr_vect, typename BlockReduce>
	static red_type red_final(red_type * rt, aggr_vect & aggr)
	{
		return BlockReduce(rt).Sum(aggr.template get<prp>());
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam aggr_vect aggregate that where i reduce
 * \tparam type of reduction to operare
 *
 */
template<typename aggr_vect, typename reduction_vectors>
struct reduce_op
{
	//! reduction aggregate
	aggr_vect & red;
	//! encapsulated input object
	reduction_vectors & input;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline reduce_op(aggr_vect & red, reduction_vectors & input)
	:red(red),input(input)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t) const
	{
		boost::mpl::at<reduction_vectors,boost::mpl::int_<T::value>>::red(red,input);
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam aggr_vect aggregate that where i reduce
 * \tparam type of reduction to operare
 *
 */
template<typename aggr_vect, typename red_type, typename reduction_vectors>
struct reduce_op_final
{
	//! reduction aggregate
	aggr_vect & red;

	//! block to reduce
	red_type * red_space;

	typedef typename boost::mpl::size<reduction_vectors>::type nred;

	red_type red_final[nred::value];

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline reduce_op_final(aggr_vect & red, red_type & red_space)
	:red(red),red_space(red_space)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		// fill the block
		red_space[threadIdx.x] = red.template get<T::value>();

		red_final[T::value] = boost::mpl::at<reduction_vectors,boost::mpl::int_<T::value>>::red_final(red_space);
	}
};


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam aggr_vect aggregate that where i reduce
 * \tparam type of reduction to operare
 *
 */
template<typename red_type, typename encap_type , typename reduction_vectors>
struct store_reduce_op_final
{
	typedef typename boost::mpl::size<reduction_vectors>::type nred;

	red_type (& red_final)[nred::value];

	encap_type & destination;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	__device__ __host__ inline store_reduce_op_final(encap_type & destination, red_type (& red_final)[nred::value])
	:red_final(red_final),destination(destination)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		typedef typename boost::mpl::at<reduction_vectors,boost::mpl::int_<T::value>>::prop prp;

		destination.template get<prp::value>() = red_final[T::value];
	}
};

template<typename vector_index_type, typename vector_index_type2, typename vector_data_type, typename reduction_operation_vector,
         typename red_type, unsigned int blockSize>
__global__ void reduce_data(vector_index_type v_offsets, vector_index_type2 v_ord, vector_data_type v_data, vector_data_type v_data_red)
{
	int offset_start = v_offsets.template get<0>(blockIdx.x);
	int offset_end = v_offsets.template get<0>(blockIdx.x+1);
	int n_block_reduce = (offset_end - offset_start) / blockDim.x;

    // Specialize BlockReduce for a 1D block
    typedef cub::BlockReduce<red_type, blockSize> BlockReduceT;

    __shared__ typename BlockReduceT::TempStorage temp_storage;

	typename vector_data_type::value_type red_v;

	typedef boost::mpl::size<reduction_operation_vector> n_reductions;

	int i = 0;
	for ( ; i < n_block_reduce ; i++)
	{
		reduce_op<typename vector_data_type::value_type,reduction_operation_vector> ro(red_v,v_data.get_o(offset_start+i*blockDim.x+threadIdx.x));

		boost::mpl::for_each_ref<boost::mpl::range_c<int,0,n_reductions::value>>(ro);
	}

	// within block reduction
	reduce_op_final<typename vector_data_type::value,red_type,reduction_operation_vector> rof(red_v,temp_storage);
	boost::mpl::for_each_ref<boost::mpl::range_c<int,0,n_reductions::value>>(rof);

	if (threadIdx.x == 0)
	{

	}
}
