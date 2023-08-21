/*
 * memory_conf.hpp
 *
 *  Created on: Aug 28, 2014
 *      Author: Pietro Incardona
 */

#ifndef MEMORY_CONF_HPP_
#define MEMORY_CONF_HPP_

#include "util/variadic_to_vmpl.hpp"
#include "t_to_memory_c.hpp"
#include "Vector/vect_isel.hpp"
#include "Vector/util.hpp"
#include "util/tokernel_transformation.hpp"
#include "util/hostDevice_util_funcs.hpp"

constexpr int SOA_layout_IA = 2;
constexpr int SOA_layout = 1;
constexpr int AOS_layout = 0;

/*! \brief This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c interleaved
 *
 * Example:
 *
 * typedef boost::mpl::vector<float,float,float[3][3], ... > A
 *
 * inter_memc<A>
 *
 * produce
 *
 * boost::fusion::vector<memory_c<float>,memory_c<float>, memory_c<multi_array<boost::mpl::vector<float,3,3>, ...... >
 *
 * \param Seq Is suppose to be an boost::mpl::vector/fusion
 *
 */
template<typename Seq>
struct inter_memc
{
	typedef typename v_transform<t_to_memory_c,Seq>::type type;
};

/*! \brief This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c<.....,MEMORY_C_REDUCED> interleaved
 *
 * This class convert a boost::mpl::fusion/vector to a boost::mpl::fusion/vector with memory_c<.....,MEMORY_C_REDUCED> interleaved
 *
 * Example:
 *
 * typedef boost::mpl::vector<float,float,float[3][3], ... > A
 *
 * inter_memc<A>
 *
 * produce
 *
 * boost::fusion::vector<memory_c<float>,memory_c<float>, memory_c<multi_array<boost::mpl::vector<float,3,3>, ...... >
 *
 * \param Seq Is suppose to be an boost::mpl::vector/fusion
 *
 */
template<typename Seq>
struct inter_memc_red
{
	typedef typename v_transform<t_to_memory_c_red,Seq>::type type;
};


/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory_traits.
 * In this implementation we interleave each property of the base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \see see inter_mem_c for detail
 *
 * \param T base type (T::type must define a boost::fusion::vector )
 *
 *
 */
template<typename T>
struct memory_traits_inte
{
	//! for each element in the vector interleave memory_c
	typedef typename inter_memc<typename T::type>::type type;

	//! indicate that it change the memory layout from the original
	typedef int yes_is_inte;

	typedef boost::mpl::int_<SOA_layout_IA> type_value;

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type, typename key_type>
	__host__ __device__ static inline auto get(data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1)))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type>
	__host__ __device__ static inline auto get_lin(data_type & data_, const g1_type & g1, size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id);
	}

	 /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a const reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type, typename key_type>
	__host__ __device__ static inline auto get_c(const data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1)))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](g1.LinId(v1));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a const reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type>
	__host__ __device__ static inline auto get_lin_c(const data_type & data_, const g1_type & g1, size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id))
	{
		return boost::fusion::at_c<p>(data_).mem_r.operator[](lin_id);
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer
	 *
	 */
	template<typename S, typename data_type, unsigned int ... prp> 
        static void hostToDevice(data_type & data_, size_t start, size_t stop)
	{
		host_to_device_impl<T,memory_traits_inte,S, prp ...> dth(data_,start,stop);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(dth);
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer
	 *
	 */
	template<typename data_type, unsigned int ... prp> 
        static void deviceToHost(data_type & data_, size_t start, size_t stop)
	{
		device_to_host_start_stop_impl<T, memory_traits_inte, prp ...> dth(data_,start,stop);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(dth);
	}
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory_traits.
 * In this implementation we interleave each property of the base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \see see inter_mem_c for detail
 *
 * \param T base type (T::type must define a boost::fusion::vector )
 *
 *
 */
/*template<typename T>
struct memory_traits_inte_red
{
	//! for each element in the vector interleave memory_c
	typedef typename inter_memc_red<typename T::type>::type type;

	//! indicate that it change the memory layout from the original
	typedef int yes_is_inte;
};*/

/*! \brief small meta-function to get the type of the memory
 *
 *
 */
template<typename T, bool is_agg>
struct memory_traits_lin_type
{
	typedef memory_c<typename T::type> type;
};

/*! \brief small meta-function to get the type of the memory
 *
 *
 */
template<typename T>
struct memory_traits_lin_type<T,false>
{
	typedef void type;
};

/*! \brief Transform the boost::fusion::vector into memory specification (memory_traits)
 *
 * Transform the boost::fusion::vector into memory specification (memory_traits).
 * In this implementation we create a buffer of base type with memory_c
 *
 * We basically create a buffer for each property
 *
 * \param T base type (T::type must define a boost::fusion::vector )
 *
 *
 */

template<typename T>
struct memory_traits_lin
{
	//! for each element in the vector interleave memory_c
	typedef typename memory_traits_lin_type<T,openfpm::vect_isel<T>::value == OPENFPM_NATIVE>::type type;

	typedef int yes_is_tlin;

	typedef boost::mpl::int_<AOS_layout> type_value;

        /*! \brief Return a reference to the selected element
         * SFINAE is used to hande CSR graphs, where boost::fusion::at_c<p>(boost::fusion::vector<>) breaks in boost 1.8* onwards
    	 *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type, typename key_type, typename std::enable_if<!std::is_same<data_type, memory_c<boost::fusion::vector<>, 1, memory>>::value,int>::type = 0>
	__host__ __device__ static inline auto get(data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)));
	}

        /*! \brief Return a reference to the selected element.
         * SFINAE is used to hande CSR graphs, where boost::fusion::at_c<p>(boost::fusion::vector<>) breaks in boost 1.8* onwards
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type, typename key_type, typename std::enable_if<std::is_same<data_type, memory_c<boost::fusion::vector<>, 1, memory>>::value,int>::type = 0>
	__host__ __device__ static inline void get(data_type & data_, const g1_type & g1, const key_type & v1) {}

        /*! \brief Return a reference to the selected element
         * SFINAE is used to hande CSR graphs, where boost::fusion::at_c<p>(boost::fusion::vector<>) breaks in boost 1.8* onwards
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
    template<unsigned int p, typename data_type, typename g1_type, typename std::enable_if<!std::is_same<data_type, memory_c<boost::fusion::vector<>, 1, memory>>::value,int>::type = 0>
	__host__ __device__ static inline auto get_lin(data_type & data_, const g1_type & g1, const size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id));
	}

        /*! \brief Return a reference to the selected element
         * SFINAE is used to hande CSR graphs, where boost::fusion::at_c<p>(boost::fusion::vector<>) breaks in boost 1.8* onwards
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
    template<unsigned int p, typename data_type, typename g1_type, typename std::enable_if<std::is_same<data_type, memory_c<boost::fusion::vector<>, 1, memory>>::value,int>::type = 0>
    __host__ __device__ static inline auto get_lin(data_type & data_, const g1_type & g1, const size_t lin_id) {}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type, typename key_type, typename std::enable_if<!std::is_same<data_type, memory_c<boost::fusion::vector<>, 1, memory>>::value,int>::type = 0>
	__host__ __device__ static inline auto get_c(const data_type & data_, const g1_type & g1, const key_type & v1) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](g1.LinId(v1)));
	}

        /*! \brief Return a reference to the selected element
         *
         * \param data object from where to take the element
         * \param g1 grid information
         * \param v1 element id
         *
         * \return a reference to the object selected
         *
         */
	template<unsigned int p, typename data_type, typename g1_type>
	__host__ __device__ static inline auto get_lin_c(const data_type & data_, const g1_type & g1, const size_t lin_id) -> decltype(boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id))) &
	{
		return boost::fusion::at_c<p>(data_.mem_r.operator[](lin_id));
	}

	/*! \brief Copy the memory from host to device
	 *
	 * \tparam (all properties are copied to prp is useless in this case)
	 *
	 * \param start start point
	 * \param stop stop point
	 *
	 */
	template<typename S, typename data_type, unsigned int ... prp> 
        static void hostToDevice(data_type & data_, size_t start, size_t stop)
	{
		data_.mem->hostToDevice(start*sizeof(T),(stop+1)*sizeof(T));
	}

	/*! \brief Synchronize the memory buffer in the device with the memory in the host
	 *
	 * \param start starting element to transfer
	 * \param stop stop element to transfer
	 *
	 * \tparam properties to transfer (ignored all properties are trasfert)
	 *
	 */
	template<typename data_type, unsigned int ... prp> 
        static void deviceToHost(data_type & data_, size_t start, size_t stop)
	{
		data_.mem->deviceToHost(start*sizeof(T),(stop+1)*sizeof(T));
	}
};


//////////////////////////////////////////////////////////////

template<typename T, typename Sfinae = void>
struct is_layout_mlin: std::false_type {};


/*! \brief is_layout_mlin
 *
 * ### Example
 *
 * \snippet util.hpp Check if the memory layout is memory_traits_lin
 *
 * return true if T is a memory_traits_lin
 *
 */
template<typename T>
struct is_layout_mlin<T, typename Void< typename T::yes_is_tlin>::type> : std::true_type
{};


template<typename T, typename Sfinae = void>
struct is_layout_inte: std::false_type {};


/*! \brief is_layout_inte
 *
 * ### Example
 *
 * \snippet util.hpp Check if the memory layout is memory_traits_inte
 *
 * return true if T is a memory_traits_inte
 *
 */
template<typename T>
struct is_layout_inte<T, typename Void< typename T::yes_is_inte>::type> : std::true_type
{};

/*! \brief is_multiple_buffer_each_prp
 *
 * return if each property is splitted on a separate buffer. This class make sense to be used if T is
 * a vector in case it is not it always return 0
 *
 * ### Example
 *
 * \snippet util_test.hpp Check if the vector has multiple buffers for each vector
 *
 *
 */
template<typename T, unsigned int impl = is_vector<T>::value >
struct is_multiple_buffer_each_prp: std::false_type
{};

template<typename T>
struct is_multiple_buffer_each_prp<T,true>: is_layout_inte<typename T::layout_base_>
{};

#endif
