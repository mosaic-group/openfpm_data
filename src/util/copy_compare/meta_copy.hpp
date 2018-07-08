#ifndef META_COPY_HPP
#define META_COPY_HPP

#include "copy_general.hpp"
#include "util/cuda_util.hpp"
#include "util/multi_array_openfpm/multi_array_ref_openfpm.hpp"

/*! \brief This class copy general objects
 *
 * * primitives
 * * array of primitives
 * * complex objects
 * * aggregates
 *
 * ### Usage of meta copy and compare for primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for primitives
 * ### Usage of meta copy and compare for array of primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for array of primitives
 * ### Usage of meta copy and compare for openfpm aggregates
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for openfpm aggregates
 * ### Usage of meta copy and compare for complex object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex object
 * ### Usage of meta copy and compare for complex aggregates object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex aggregates object
 * ### Usage of meta copy and compare for Point_test
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for Point_test
 *
 */
template<typename T>
struct meta_copy
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	__device__ __host__ static inline void meta_copy_(const T & src, T & dst)
	{
		copy_general<T>(src,dst);
	}

	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	__device__ __host__  static inline void meta_copy_(const T & src, T && dst)
	{
		copy_general<T>(src,dst);
	}
};

/*! \brief copy for a source object to a destination
 *
 * \tparam Tsrc source object
 * \tparam Tdst destination object
 *
 */
template<typename Tsrc,typename Tdst>
struct meta_copy_d
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_d_(const Tsrc & src, Tdst & dst)
	{
		copy_general<Tsrc>(src,dst);
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct meta_copy<T[N1]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	__device__ __host__ static inline void meta_copy_(const T src[N1], T dst[N1])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			copy_general<T>(src[i1],dst[i1]);
		}
	}
};

//! Partial specialization for N=1 1D-Array
template<typename Tsrc, typename Tdst, size_t N1>
struct meta_copy_d<Tsrc[N1],Tdst>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_d_(const Tsrc src[N1], Tdst && dst)
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			copy_general<Tsrc>(src[i1],static_cast<Tsrc&>(dst[i1]));
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct meta_copy<T[N1][N2]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	__device__ __host__ static inline void meta_copy_(const T src[N1][N2], T dst[N1][N2])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				copy_general<T>(src[i1][i2],dst[i1][i2]);
			}
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename Tsrc, typename Tdst,size_t N1,size_t N2>
struct meta_copy_d<Tsrc[N1][N2],Tdst>
{
	static inline void meta_copy_d_(const Tsrc src[N1][N2], Tdst && dst)
	{
		/*! \brief copy and object from src to dst
		 *
		 * \param src source object to copy
		 * \param dst destination object
		 *
		 */
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				copy_general<Tsrc>(src[i1][i2],static_cast<Tsrc&>(dst[i1][i2]));
			}
		}
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct meta_copy<T[N1][N2][N3]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3], T dst[N1][N2][N3])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					copy_general<T>(src[i1][i2][i3],dst[i1][i2][i3]);
				}
			}
		}
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct meta_copy<T[N1][N2][N3][N4]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4], T dst[N1][N2][N3][N4])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						copy_general<T>(src[i1][i2][i3][i4],dst[i1][i2][i3][i4]);
					}
				}
			}
		}
	}
};

//! Partial specialization for N=5
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5>
struct meta_copy<T[N1][N2][N3][N4][N5]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5], T dst[N1][N2][N3][N4][N5])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							copy_general<T>(src[i1][i2][i3][i4][i5],dst[i1][i2][i3][i4][i5]);
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=6
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6>
struct meta_copy<T[N1][N2][N3][N4][N5][N6]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5][N6], T dst[N1][N2][N3][N4][N5][N6])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								copy_general<T>(src[i1][i2][i3][i4][i5][i6],dst[i1][i2][i3][i4][i5][i6]);
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=7
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5][N6][N7], T dst[N1][N2][N3][N4][N5][N6][N7])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									copy_general<T>(src[i1][i2][i3][i4][i5][i6][i7],dst[i1][i2][i3][i4][i5][i6][i7]);
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=8
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5][N6][N7][N8], T dst[N1][N2][N3][N4][N5][N6][N7][N8])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										copy_general<T>(src[i1][i2][i3][i4][i5][i6][i7][i8],dst[i1][i2][i3][i4][i5][i6][i7][i8]);
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=9
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8][N9]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5][N6][N7][N8][N9], T dst[N1][N2][N3][N4][N5][N6][N7][N8][N9])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										for (size_t i9 = 0 ; i9 < N9 ; i9++)
										{
											copy_general<T>(src[i1][i2][i3][i4][i5][i6][i7][i8][i9],dst[i1][i2][i3][i4][i5][i6][i7][i8][i9]);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=10
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9, size_t N10>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10]>
{
	/*! \brief copy and object from src to dst
	 *
	 * \param src source object to copy
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_(const T src[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10], T dst[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										for (size_t i9 = 0 ; i9 < N9 ; i9++)
										{
											for (size_t i10 = 0 ; i10 < N10 ; i10++)
											{
												copy_general<T>(src[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10],dst[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10]);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

////////////////////////// VERSION WITH OPERATION ////////////////////////


/*! \brief This class copy general objects applying an operation
 *
 * * primitives
 * * array of primitives
 * * complex objects
 * * aggregates
 *
 * ### Usage of meta copy and compare for primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for primitives
 * ### Usage of meta copy and compare for array of primitives
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for array of primitives
 * ### Usage of meta copy and compare for openfpm aggregates
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for openfpm aggregates
 * ### Usage of meta copy and compare for complex object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex object
 * ### Usage of meta copy and compare for complex aggregates object
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for complex aggregates object
 * ### Usage of meta copy and compare for Point_test
 * \snippet meta_cc_unit_tests.hpp Usage of meta copy and compare for Point_test
 *
 */
template<template<typename,typename> class op, typename T>
struct meta_copy_op
{
	/*! \brief Meta-copy applying an operation
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_op_(const T & src, T & dst)
	{
		copy_general_op<op,T>(src,dst);
	}

	/*! \brief Meta-copy applying an operation
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_op_(const T & src, T && dst)
	{
		copy_general_op<op,T>(src,dst);
	}
};


//! Partial specialization for N=1 1D-Array
template<template<typename,typename> class op, typename T,size_t N1>
struct meta_copy_op<op,T[N1]>
{
	/*! \brief Meta-copy applying an operation
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_op_(const T src[N1], T dst[N1])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			copy_general_op<op,T>(src[i1],dst[i1]);
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<template<typename,typename> class op, typename T,size_t N1,size_t N2>
struct meta_copy_op<op,T[N1][N2]>
{
	/*! \brief Meta-copy applying an operation
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_op_(const T src[N1][N2], T dst[N1][N2])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				copy_general_op<op,T>(src[i1][i2],dst[i1][i2]);
			}
		}
	}
};


//! Partial specialization for N=3
template<template<typename,typename> class op, typename T,size_t N1,size_t N2,size_t N3>
struct meta_copy_op<op,T[N1][N2][N3]>
{
	/*! \brief Meta-copy applying an operation
	 *
	 * \param src source object
	 * \param dst destination object
	 *
	 */
	static inline void meta_copy_op_(const T src[N1][N2][N3], T dst[N1][N2][N3])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					copy_general_op<op,T>(src[i1][i2][i3],dst[i1][i2][i3]);
				}
			}
		}
	}
};

#endif
