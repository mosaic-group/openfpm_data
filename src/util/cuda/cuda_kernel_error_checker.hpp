/*
 * se_class1_cuda.hpp
 *
 *  Created on: Jan 13, 2019
 *      Author: i-bird
 */

#ifndef SE_CLASS1_CUDA_HPP_
#define SE_CLASS1_CUDA_HPP_

#include "Grid/util.hpp"
#include "Vector/util.hpp"

template<typename T, int type_of_t=2*is_grid<T>::value+1*is_vector<T>::value>
struct check_type
{
	static int check(void * ptr, int prp, T & arg)
	{
		return false;
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. It check if the
 * pointer ptr match one of the pointer properties
 *
 */
template<typename data_type>
struct check_device_ptr
{
	//! pointer to check
	void * ptr;

	//! Data to check
	data_type & data;

	int prp;

	mutable bool result;

	/*! \brief constructor
	 *
	 * \param ptr pointer to check
	 * \param data data structure
	 *
	 */
	inline check_device_ptr(void * ptr, int prp, data_type & data)
	:ptr(ptr),data(data),prp(prp)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		if (T::value == prp)
		{
			result = data.template getPointer<T::value>() == ptr;
		}
	}
};

template<typename T>
struct check_type<T,1>
{
	static int check(void * ptr, int prp, T & arg)
	{
		check_device_ptr<T> cp(ptr,prp,arg);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::value_type::max_prop> >(cp);

		return cp.result;
	}
};

template<typename T>
struct check_type<T,2>
{
	static int check(void * ptr, int prp, T & arg)
	{
		check_device_ptr<T> cp(ptr,prp,arg);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::value_type::max_prop> >(cp);

		return cp.result;
	}
};

template<typename ArgL>
int error_args_impl(void * ptr, int prp, ArgL argl)
{
	if (check_type<ArgL>::check(ptr,prp,argl) == true)
	{
		return 0;
	}
	return -1;
}

template<typename ArgL, typename ... Args>
int error_args_impl(void * ptr, int prp, ArgL argl, Args ... args)
{
	if (check_type<ArgL>::check(ptr,prp,argl) == true)
	{
		return sizeof...(args);
	}
	return error_args_impl(ptr, prp, args ...);
}

template<typename ... Args>int error_arg(void * ptr, int prp, Args ... args)
{
	int pos = error_args_impl(ptr, prp, args ... );
	return sizeof...(args) - pos - 1;
}

#include <boost/algorithm/string.hpp>

#ifdef SE_CLASS1
#define CHECK_SE_CLASS1_PRE int dev_mem[] = {0,0,0,0,0,0,0,0,0,0,0};
//#define CHECK_SE_CLASS1_POST(...) cudaMemcpyFromSymbol(dev_mem,global_cuda_error_array,sizeof(dev_mem));
#define CHECK_SE_CLASS1_POST(kernel_call,...) cudaMemcpyFromSymbol(dev_mem,global_cuda_error_array,sizeof(dev_mem)); \
		                     if (dev_mem[0] != 0)\
		                     {\
		                    	 void * ptr = (void *)*(size_t *)&dev_mem[1]; \
		                    	 int prp_err = dev_mem[3];\
		                    	 int ea = error_arg(ptr,prp_err,__VA_ARGS__);\
		                    	 std::string args_s( #__VA_ARGS__ );\
		                    	 std::vector<std::string> results;\
		                    	 boost::split(results, args_s, [](char c){return c == ',';});\
		                    	 std::cout << __FILE__ << ":" << __LINE__ << " Overflow detected in Kernel: " << kernel_call << " from the structure " << results[ea] << " property: " << prp_err << " index:(" ;\
		                    	 int i = 0; \
		                    	 for ( ; i < dev_mem[4]-1 ; i++)\
		                    	 {\
		                    		 std::cout << dev_mem[5+i] << ",";\
								 }\
								 std::cout << dev_mem[5+i];\
								 std::cout << ")";\
								 std::cout << " thread: " << "(" << dev_mem[6+i] << "," << dev_mem[7+i] << "," << dev_mem[8+i] << ")*(" << dev_mem[9+i] << "," << dev_mem[10+i] << "," << dev_mem[11+i] << ")+(" << dev_mem[12+i] << "," << dev_mem[13+i] << "," << dev_mem[14+i] << ")" << std::endl;\
		                     }
#else
#define CHECK_SE_CLASS1_PRE
#define CHECK_SE_CLASS1_POST(kernel_call,...)
#endif


#endif /* SE_CLASS1_CUDA_HPP_ */
