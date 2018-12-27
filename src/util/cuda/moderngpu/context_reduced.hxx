/*
 * context_reduced.hxx
 *
 *  Created on: Dec 27, 2018
 *      Author: i-bird
 */

#ifndef CONTEXT_REDUCED_HXX_
#define CONTEXT_REDUCED_HXX_

#include <cstdarg>
#include <string>

namespace mgpu {

enum memory_space_t {
  memory_space_device = 0,
  memory_space_host = 1
};

struct cuda_exception_t : std::exception {
  cudaError_t result;

  cuda_exception_t(cudaError_t result_) : result(result_) { }
  virtual const char* what() const noexcept {
    return cudaGetErrorString(result);
  }
};

namespace detail {

inline std::string stringprintf(const char* format, ...) {
  va_list args;
  va_start(args, format);
  int len = vsnprintf(0, 0, format, args);
  va_end(args);

  // allocate space.
  std::string text;
  text.resize(len);

  va_start(args, format);
  vsnprintf(&text[0], len + 1, format, args);
  va_end(args);

  return text;
}

} // namespace detail

inline std::string device_prop_string(cudaDeviceProp prop) {
  int ordinal;
  cudaGetDevice(&ordinal);

  size_t freeMem, totalMem;
  cudaError_t result = cudaMemGetInfo(&freeMem, &totalMem);
  if(cudaSuccess != result) throw cuda_exception_t(result);

  double memBandwidth = (prop.memoryClockRate * 1000.0) *
    (prop.memoryBusWidth / 8 * 2) / 1.0e9;

  std::string s = detail::stringprintf(
    "%s : %8.3lf Mhz   (Ordinal %d)\n"
    "%d SMs enabled. Compute Capability sm_%d%d\n"
    "FreeMem: %6dMB   TotalMem: %6dMB   %2d-bit pointers.\n"
    "Mem Clock: %8.3lf Mhz x %d bits   (%5.1lf GB/s)\n"
    "ECC %s\n\n",
    prop.name, prop.clockRate / 1000.0, ordinal,
    prop.multiProcessorCount, prop.major, prop.minor,
    (int)(freeMem / (1<< 20)), (int)(totalMem / (1<< 20)), 8 * sizeof(int*),
    prop.memoryClockRate / 1000.0, prop.memoryBusWidth, memBandwidth,
    prop.ECCEnabled ? "Enabled" : "Disabled");
  return s;
}

////////////////////////////////////////////////////////////////////////////////
// context_t
// Derive context_t to add support for streams and a custom allocator.

struct context_t {
  context_t() = default;

  // Disable copy ctor and assignment operator. We don't want to let the
  // user copy only a slice.
  context_t(const context_t& rhs) = delete;
  context_t& operator=(const context_t& rhs) = delete;

  virtual const cudaDeviceProp& props() const = 0;
  virtual int ptx_version() const = 0;
  virtual cudaStream_t stream() = 0;

  // Alloc GPU memory.
  virtual void* alloc(size_t size, memory_space_t space) = 0;
  virtual void free(void* p, memory_space_t space) = 0;

  // cudaStreamSynchronize or cudaDeviceSynchronize for stream 0.
  virtual void synchronize() = 0;

  virtual cudaEvent_t event() = 0;
  virtual void timer_begin() = 0;
  virtual double timer_end() = 0;
};

}

#endif /* CONTEXT_REDUCED_HXX_ */
