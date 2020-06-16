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

struct pointer_check
{
	//! Indicate if the pointer match
	bool match;

	//! match string
	std::string match_str;
};

template<typename T, int type_of_t=has_check_device_pointer<T>::value>
struct check_type
{
	static pointer_check check(void * ptr, int prp, T & arg)
	{
		pointer_check pc;

		pc.match = false;

		return pc;
	}
};



template<typename T>
struct check_type<T,1>
{
	static pointer_check check(void * ptr, int prp, T & arg)
	{
		return arg.check_device_pointer(ptr);
	}
};

struct pos_pc
{
	int pos;
	pointer_check pc;
};

template<typename ArgL>
pos_pc error_args_impl(void * ptr, int prp, ArgL argl)
{
	pos_pc pp;
	pointer_check pc = check_type<ArgL>::check(ptr,prp,argl);
	if (pc.match == true)
	{
		pp.pos = 0;
		pp.pc = pc;
		return pp;
	}

	pp.pos = -1;

	return pp;
}

template<typename ArgL, typename ... Args>
pos_pc error_args_impl(void * ptr, int prp, ArgL argl, Args ... args)
{
	pos_pc pp;
	pointer_check pc = check_type<ArgL>::check(ptr,prp,argl);
	if (pc.match == true)
	{
		pp.pos = sizeof...(args);
		pp.pc = pc;
		return pp;
	}
	return error_args_impl(ptr, prp, args ...);
}

template<typename ... Args>pos_pc error_arg(void * ptr, int prp, Args ... args)
{
	pos_pc pp;
	pp = error_args_impl(ptr, prp, args ... );
	pp.pos = sizeof...(args) - pp.pos - 1;
	return pp;
}

#include <boost/algorithm/string.hpp>

#ifdef SE_CLASS1
#define CUDA_LAUNCH_ERROR_OBJECT std::runtime_error("Runtime vector error");
#define CHECK_SE_CLASS1_PRE int dev_mem[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#define CHECK_SE_CLASS1_POST(kernel_call,...) cudaMemcpyFromSymbol(dev_mem,global_cuda_error_array,sizeof(dev_mem)); \
		                     if (dev_mem[0] != 0)\
		                     {\
		                    	 void * ptr = (void *)*(size_t *)&dev_mem[1]; \
		                    	 int prp_err = dev_mem[3];\
		                    	 pos_pc ea = error_arg(ptr,prp_err,__VA_ARGS__);\
		                    	 std::string args_s( #__VA_ARGS__ );\
		                    	 std::vector<std::string> results;\
		                    	 boost::split(results, args_s, [](char c){return c == ',';});\
		                    	 std::string data_s;\
		                    	 if (ea.pos >= results.size())\
		                    	 {data_s = "Internal";}\
								 else\
								 {data_s = results[ea.pos];}\
		                    	 std::cout << __FILE__ << ":" << __LINE__ << " Overflow detected in Kernel: " << kernel_call << " from the structure: " << data_s << " property: " << prp_err << " index:(" ;\
		                    	 int i = 0; \
		                    	 for ( ; i < dev_mem[4]-1 ; i++)\
		                    	 {\
		                    		 std::cout << dev_mem[5+i] << ",";\
								 }\
								 std::cout << dev_mem[5+i];\
								 std::cout << ")";\
								 std::cout << " thread: " << "(" << dev_mem[6+i] << "," << dev_mem[7+i] << "," << dev_mem[8+i] << ")*(" << dev_mem[9+i] << "," << dev_mem[10+i] << "," << dev_mem[11+i] << ")+(" << dev_mem[12+i] << "," << dev_mem[13+i] << "," << dev_mem[14+i] << ")" << std::endl;\
		                    	 std::cout << "Internal error report: " << ea.pc.match_str << std::endl;\
								 ACTION_ON_ERROR(CUDA_LAUNCH_ERROR_OBJECT);\
		                     }
#else
#define CHECK_SE_CLASS1_PRE
#define CHECK_SE_CLASS1_POST(kernel_call,...)
#endif


#endif /* SE_CLASS1_CUDA_HPP_ */
