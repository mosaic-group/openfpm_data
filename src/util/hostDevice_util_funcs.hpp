#ifndef HOSTDEVICE_UTIL_FUNCS_HPP_
#define HOSTDEVICE_UTIL_FUNCS_HPP_

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename T_type, template<typename> class layout_base, unsigned int ... prp>
struct device_to_host_start_stop_impl
{
	//! encapsulated destination object
	typename layout_base<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! start
	size_t start;

	//! stop
	size_t stop;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline device_to_host_start_stop_impl(typename layout_base<T_type>::type & dst,size_t start,size_t stop)
	:dst(dst),start(start),stop(stop)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type prp_id;

		typedef typename boost::mpl::at<typename T_type::type,prp_id>::type p_type;

		boost::fusion::at_c<prp_id::value>(dst).mem->deviceToHost(start*sizeof(p_type),(stop+1)*sizeof(p_type));
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename T_type, template<typename> class layout_base, unsigned int ... prp>
struct device_to_host_impl
{
	//! encapsulated destination object
	typename layout_base<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline device_to_host_impl(typename layout_base<T_type>::type & dst)
	:dst(dst)
	{
	};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem->deviceToHost();
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to copy one encap into another encap object
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename T_type, template<typename> class layout_base , typename Memory, unsigned int ... prp>
struct host_to_device_impl
{
	//! encapsulated destination object
	typename layout_base<T_type>::type & dst;

	//! Convert the packed properties into an MPL vector
	typedef typename to_boost_vmpl<prp...>::type v_prp;

	//! starting element
	size_t start;

	//! stop element
	size_t stop;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline host_to_device_impl(typename layout_base<T_type>::type & dst,size_t start, size_t stop)
	:dst(dst),start(start),stop(stop)
	{};


	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef typename boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type ele_type;

		typedef decltype(boost::fusion::at_c<ele_type::value>(dst).mem_r) mem_r_type;

		typedef typename boost::mpl::at<typename T_type::type,ele_type>::type type_prp;

		typedef typename toKernel_transform<layout_base,typename mem_r_type::value_type>::type kernel_type;

		typedef boost::mpl::int_<(is_vector<typename mem_r_type::value_type>::value ||
								  is_vector_dist<typename mem_r_type::value_type>::value ||
								  is_gpu_celllist<typename mem_r_type::value_type>::value) + 2*std::is_array<type_prp>::value + std::rank<type_prp>::value> crh_cond;

		call_recursive_host_device_if_vector<typename mem_r_type::value_type,
											 kernel_type,
											 type_prp,
											 layout_base,
											 crh_cond::value>
		::template transform<Memory,mem_r_type>(static_cast<Memory *>(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem),
									 boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem_r,
				                       start*sizeof(type_prp),
				                       (stop+1)*sizeof(type_prp));

		// here we have to recursively call hostToDevice for each nested vector
		call_recursive_host_device_if_vector<typename mem_r_type::value_type,
											 kernel_type,
											 type_prp,
											 layout_base,
											 0>
		::call(boost::fusion::at_c<boost::mpl::at<v_prp,boost::mpl::int_<T::value>>::type::value>(dst).mem_r,start,stop);
	}
};

#endif