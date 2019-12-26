/*
 * ofp_context.hxx
 *
 *  Created on: Nov 15, 2018
 *      Author: i-bird
 */

#ifndef OFP_CONTEXT_HXX_
#define OFP_CONTEXT_HXX_

#include <iostream>

#ifdef CUDA_GPU

	#if defined(__NVCC__) || defined(__HIPCC__)
	#include "util/cuda/moderngpu/context.hxx"

	namespace mgpu
	{
		enum gpu_context_opt
		{
			no_print_props,//!< no_print_props
			print_props,   //!< print_props
			dummy          //!< dummy
		};


		////////////////////////////////////////////////////////////////////////////////
		// standard_context_t is a trivial implementation of context_t. Users can
		// derive this type to provide a custom allocator.

		class ofp_context_t : public context_t
		{
			protected:
				cudaDeviceProp _props;
				int _ptx_version;
				cudaStream_t _stream;

				cudaEvent_t _timer[2];
				cudaEvent_t _event;

				openfpm::vector_gpu<aggregate<unsigned char>> tmem;

				// Making this a template argument means we won't generate an instance
				// of dummy_k for each translation unit.
				template<int dummy_arg = 0>
				void init(int dev_num, gpu_context_opt opt)
				{
					cudaFuncAttributes attr;
					cudaError_t result = cudaFuncGetAttributes(&attr, dummy_k<0>);
					if(cudaSuccess != result) throw cuda_exception_t(result);
					_ptx_version = attr.ptxVersion;

					int num_dev;
					cudaGetDeviceCount(&num_dev);

					if (num_dev == 0) {return;}

					if (opt != gpu_context_opt::dummy)
					{
						cudaSetDevice(dev_num % num_dev);
					}

					int ord;
					cudaGetDevice(&ord);
					cudaGetDeviceProperties(&_props, ord);

					cudaEventCreate(&_timer[0]);
					cudaEventCreate(&_timer[1]);
					cudaEventCreate(&_event);
				}

			public:


				/*! \brief gpu context constructor
				 *
				 * \param opt options for this gpu context
				 *
				 */
				ofp_context_t(gpu_context_opt opt = gpu_context_opt::no_print_props , int dev_num = 0, cudaStream_t stream_ = 0)
				:context_t(), _stream(stream_)
				{
					init(dev_num,opt);
					if(opt == gpu_context_opt::print_props)
					{
						printf("%s\n", device_prop_string(_props).c_str());
					}
				}

				~ofp_context_t()
				{
					cudaEventDestroy(_timer[0]);
					cudaEventDestroy(_timer[1]);
					cudaEventDestroy(_event);
				}

				virtual const cudaDeviceProp& props() const { return _props; }
				virtual int ptx_version() const { return _ptx_version; }
				virtual cudaStream_t stream() { return _stream; }

				// Alloc GPU memory.
				virtual void* alloc(size_t size, memory_space_t space)
				{
					void* p = nullptr;
					if(size)
					{
						cudaError_t result = (memory_space_device == space) ?cudaMalloc(&p, size) : cudaMallocHost(&p, size);
						if(cudaSuccess != result) throw cuda_exception_t(result);
					}
					return p;
				}

				virtual void free(void* p, memory_space_t space)
				{
					if(p)
					{
						cudaError_t result = (memory_space_device == space) ? cudaFree(p) : cudaFreeHost(p);
						if(cudaSuccess != result) throw cuda_exception_t(result);
					}
				}

				virtual void synchronize()
				{
					cudaError_t result = _stream ?
					cudaStreamSynchronize(_stream) :
					cudaDeviceSynchronize();
					if(cudaSuccess != result) throw cuda_exception_t(result);
				}

				virtual cudaEvent_t event()
				{
					return _event;
				}

				virtual void timer_begin()
				{
					cudaEventRecord(_timer[0], _stream);
				}

				virtual double timer_end()
				{
					cudaEventRecord(_timer[1], _stream);
					cudaEventSynchronize(_timer[1]);
					float ms;
					cudaEventElapsedTime(&ms, _timer[0], _timer[1]);
					return ms / 1.0e3;
				}

				virtual int getDevice()
				{
					int dev = 0;

					cudaGetDevice(&dev);

					return dev;
				}

				virtual int getNDevice()
				{
					int num_dev;
					cudaGetDeviceCount(&num_dev);

					return num_dev;
				}

				openfpm::vector_gpu<aggregate<unsigned char>> & getTemporalCUB()
				{
					return tmem;
				}
		};

	}

	#else

	#include "util/cuda/moderngpu/context_reduced.hxx"

	namespace mgpu
	{
		enum gpu_context_opt
		{
			no_print_props,//!< no_print_props
			print_props,   //!< print_props
			dummy          //!< dummy
		};


		////////////////////////////////////////////////////////////////////////////////
		// standard_context_t is a trivial implementation of context_t. Users can
		// derive this type to provide a custom allocator.

		class ofp_context_t : public context_t
		{
			protected:
				cudaDeviceProp _props;
				int _ptx_version;
				cudaStream_t _stream;

				cudaEvent_t _timer[2];
				cudaEvent_t _event;

				// Making this a template argument means we won't generate an instance
				// of dummy_k for each translation unit.
				template<int dummy_arg = 0>
				void init(int dev_num, gpu_context_opt opt)
				{
					cudaFuncAttributes attr;

					_ptx_version = 0;

					int num_dev;
					cudaGetDeviceCount(&num_dev);

					if (num_dev == 0) {return;}

					if (opt != gpu_context_opt::dummy)
					{
						cudaSetDevice(dev_num % num_dev);
					}

					int ord;
					cudaGetDevice(&ord);
					cudaGetDeviceProperties(&_props, ord);

					cudaEventCreate(&_timer[0]);
					cudaEventCreate(&_timer[1]);
					cudaEventCreate(&_event);
				}

			public:

				/*! \brief gpu context constructor
				 *
				 * \param opt options for this gpu context
				 *
				 */
				ofp_context_t(gpu_context_opt opt = gpu_context_opt::no_print_props , int dev_num = 0, cudaStream_t stream_ = 0)
				:context_t(), _stream(stream_)
				{
					init(dev_num,opt);
					if(opt == gpu_context_opt::print_props)
					{
						printf("%s\n", device_prop_string(_props).c_str());
					}
				}

				~ofp_context_t()
				{
					cudaEventDestroy(_timer[0]);
					cudaEventDestroy(_timer[1]);
					cudaEventDestroy(_event);
				}

				virtual const cudaDeviceProp& props() const
				{
					return _props;
				}

				virtual int ptx_version() const
				{
					std::cout << __FILE__ << ":" << __LINE__ << " error to use this function you must compile the class ofp_context_t with NVCC" << std::endl;
					return 0;
				}

				virtual cudaStream_t stream() { return _stream; }

				// Alloc GPU memory.
				virtual void* alloc(size_t size, memory_space_t space)
				{
					void* p = nullptr;
					if(size)
					{
						cudaError_t result = (memory_space_device == space) ?cudaMalloc(&p, size) : cudaMallocHost(&p, size);
						if(cudaSuccess != result) throw cuda_exception_t(result);
					}
					return p;
				}

				virtual void free(void* p, memory_space_t space)
				{
					if(p)
					{
						cudaError_t result = (memory_space_device == space) ? cudaFree(p) : cudaFreeHost(p);
						if(cudaSuccess != result) throw cuda_exception_t(result);
					}
				}

				virtual void synchronize()
				{
					cudaError_t result = _stream ?
					cudaStreamSynchronize(_stream) :
					cudaDeviceSynchronize();
					if(cudaSuccess != result) throw cuda_exception_t(result);
				}

				virtual cudaEvent_t event()
				{
					return _event;
				}

				virtual void timer_begin()
				{
					cudaEventRecord(_timer[0], _stream);
				}

				virtual double timer_end()
				{
					cudaEventRecord(_timer[1], _stream);
					cudaEventSynchronize(_timer[1]);
					float ms;
					cudaEventElapsedTime(&ms, _timer[0], _timer[1]);
					return ms / 1.0e3;
				}

				virtual int getDevice()
				{
					int dev = 0;

					cudaGetDevice(&dev);

					return dev;
				}
		};

	}

	#endif

#else

	namespace mgpu
	{

		enum gpu_context_opt
		{
			no_print_props,//!< no_print_props
			print_props,   //!< print_props
			dummy          //!< dummy
		};

		// Stub class for modern gpu

		struct ofp_context_t
		{
			ofp_context_t(gpu_context_opt opt = gpu_context_opt::no_print_props , int dev_num = 0)
			{}
		};
	}

#endif


#endif /* OFP_CONTEXT_HXX_ */
