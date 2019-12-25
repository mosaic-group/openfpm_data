//
// Created by tommaso on 24/05/19.
//

#ifndef OPENFPM_PDATA_BLOCKMAPGPU_DIMENSIONALITYWRAPPERS_CUH
#define OPENFPM_PDATA_BLOCKMAPGPU_DIMENSIONALITYWRAPPERS_CUH

#include <type_traits>

// Type management and wrapping

/**
 * Get the type of the array of scalars from the block array type
 *
 * @tparam BaseT The block array type
 */
template<typename BaseT>
struct ComposeArrayType
{
    typedef typename BaseT::scalarType type;
};

template<typename BaseType, unsigned int N1>
struct ComposeArrayType<BaseType[N1]>
{
    typedef typename BaseType::scalarType type[N1];
};

template<typename BaseType, unsigned int N1, unsigned int N2>
struct ComposeArrayType<BaseType[N1][N2]>
{
    typedef typename BaseType::scalarType type[N1][N2];
};

// Data management and wrapping

/**
 * Provides a view on data from an array of blocks,
 * yielding items at the same offset within blocks through the block array.
 *
 * @tparam BaseT
 * @tparam Nup
 * @tparam N
 */
template<typename BaseT, unsigned int Nup, unsigned int N>
class MultiArrayViewGpu
{
    typedef typename std::remove_extent<BaseT>::type array_slice;
    typedef std::extent<BaseT> ext;

    BaseT * ptr;

public:
    __device__ __host__ inline MultiArrayViewGpu(BaseT * ptr)
            :ptr(ptr)
    {}

    __device__ __host__ inline MultiArrayViewGpu<array_slice, ext::value, N-1> operator[](int i)
    {
        return MultiArrayViewGpu< array_slice,ext::value,N-1>((array_slice *)(ptr+i));
    }
};

template<typename BaseT, unsigned int Nup>
class MultiArrayViewGpu<BaseT,Nup,0>
{
    BaseT * ptr;

public:

    __device__ __host__ MultiArrayViewGpu(BaseT * ptr)
            :ptr(ptr)
    {}

    template <typename IndexT>
    __device__ __host__ MultiArrayViewGpu(BaseT * ptr, IndexT offset)
            :ptr((BaseT*)(((typename BaseT::scalarType *)ptr) + offset))
    {}

    __device__ __host__ inline typename BaseT::scalarType & operator[](int i)
    {
        return *((typename BaseT::scalarType *)(&(ptr[i])));
    }

    __device__ __host__ inline typename BaseT::scalarType & operator[](int i) const
    {
        return *((typename BaseT::scalarType *)(&(ptr[i])));
    }

    template <typename T>
    __device__ __host__ inline MultiArrayViewGpu<BaseT,Nup,0> & operator=(const T& other)
    {
        for (int i=0; i< Nup ; i++)
        {
            this->operator[](i) = other[i];
        }
        return *this;
    }
};

/**
 * Provides an on-stack type wrapping the original type and switching index access
 * @tparam BaseT
 */
template<typename BaseT>
struct ArrayWrapper
{
    BaseT data;

    __device__ __host__ inline typename BaseT::scalarType & operator[](int i)
    {
        return data[i];
    }

    __device__ __host__ inline const typename BaseT::scalarType & operator[](int i) const
    {
        return data[i];
    }
};

template<typename BaseType, unsigned int N1>
struct ArrayWrapper<BaseType[N1]>
{
    BaseType data[N1];

    __device__ __host__ MultiArrayViewGpu<BaseType,N1,0> operator[](int i)
    {
        return MultiArrayViewGpu<BaseType,N1,0>((BaseType*)(((typename BaseType::scalarType *)data) + i));
    }

    __device__ __host__ const MultiArrayViewGpu<BaseType,N1,0> operator[](int i) const
    {
        return MultiArrayViewGpu<BaseType,N1,0>((BaseType*)(((typename BaseType::scalarType *)data) + i));
    }
};

template<typename BaseType, unsigned int N1, unsigned int N2>
struct ArrayWrapper<BaseType[N1][N2]>
{
    BaseType array[N1][N2];

    __device__ __host__ MultiArrayViewGpu<BaseType[N2],N1,1> operator[](int i)
    {
        return MultiArrayViewGpu<BaseType[N2],N1,1>((BaseType(*)[N2])(((typename BaseType::scalarType *)array) + i));
    }
};

/**
 * A wrapper for general dimension right hand sides for assignments.
 * This involves copies of data into local registers upon construction.
 *
 * @tparam BlockT
 */
template<typename BlockT>
struct RhsBlockWrapper
{
    typename BlockT::scalarType &value;

    template <typename IndexT>
    __device__ __host__ RhsBlockWrapper(BlockT &block, IndexT offset) : value(block[offset]) {}
};

template<typename BlockT, unsigned int N>
struct RhsBlockWrapper<BlockT[N]>
{
    typename BlockT::scalarType value[N];

    template <typename T, typename IndexT>
    __device__ __host__ RhsBlockWrapper(T input, IndexT offset)
    {
        for (int i=0; i<N; ++i)
        {
            value[i] = input[i][offset];
        }
    }
};

template<typename BlockT, unsigned int N1, unsigned int N2>
struct RhsBlockWrapper<BlockT[N1][N2]>
{
    typename BlockT::scalarType value[N1][N2];

    template <typename T, typename IndexT>
    __device__ __host__ RhsBlockWrapper(T input, IndexT offset)
    {
        for (int i=0; i<N1; ++i)
        {
            for (int j=0; j<N2; ++j)
            {
                value[i][j] = input[i][j][offset];
            }
        }
    }
};

/**
 * Functor which provides primitives to be used to deal with general dimensional objects
 *
 * @tparam T1 The type providing dimensionality information
 */
template <typename T1>
struct generalDimensionFunctor
{
    template <typename T1b, typename T2, typename T3>
    __device__ __host__ inline static void assignWithOffsetRHS(T1b &dst, const T2 &src, T3 offset)
    {
        dst = src[offset];
    }

    template <typename T1b, typename T2, typename T3>
    __device__ __host__ inline static void assignWithOffset(T1b &dst, const T2 &src, T3 offset)
    {
        dst[offset] = src[offset];
    }

    template<typename op, typename T1b, typename T2>
    __device__ inline static void applyOp(T1b &a, const T2 &b, bool aExist, bool bExist)
    {
        op op_;
        if (aExist && bExist)
        {
            a = op_(a, b);
        }
        else if (bExist)
        {
            a = b;
        }
    }
};

template <typename T1, unsigned int N1>
struct generalDimensionFunctor<T1[N1]>
{
    template <typename T1b, typename T2, typename T3>
    __device__ __host__ inline static void assignWithOffsetRHS(T1b (& dst)[N1], const T2 &src, T3 offset)
    {
        for (int i = 0; i < N1; ++i)
        {
//            dst[i] = src[i][offset];
            generalDimensionFunctor<T1>::assignWithOffsetRHS(dst[i], src[i], offset);
        }
    }

    template <typename T1b, typename T2, typename T3>
//    __device__ __host__ inline static void assignWithOffset(T1b (& dst)[N1], const T2 &src, T3 offset)
    __device__ __host__ inline static void assignWithOffset(T1b dst, const T2 &src, T3 offset)
    {
        for (int i = 0; i < N1; ++i)
        {
//            dst[i][offset] = src[i][offset];
            generalDimensionFunctor<T1>::assignWithOffset(dst[i], src[i], offset);
        }
    }

    template<typename op, typename T1b, typename T2>
    __device__ inline static void applyOp(T1b &a, const T2 &b, bool aExist, bool bExist)
    {
        for (int i = 0; i < N1; ++i)
        {
            generalDimensionFunctor<T1>::template applyOp<op>(a[i], b[i], aExist, bExist);
        }
    }

    template<typename op, typename T1b, typename T2>
    __device__ inline static void applyOp(const T1b &a, const T2 &b, bool aExist, bool bExist)
    {
        for (int i = 0; i < N1; ++i)
        {
            generalDimensionFunctor<T1>::template applyOp<op>(a[i], b[i], aExist, bExist);
        }
    }
};

//template <typename T1, unsigned int N1, unsigned int N2>
//struct generalDimensionFunctor<T1[N1][N2]>
//{
//    template <typename T2, typename T3>
//    __device__ __host__ inline static void assign(T1 (& dst)[N1][N2], const T2 &src, T3 offset)
//    {
//        for (int i = 0; i < N1; ++i)
//        {
//            for (int j = 0; j < N2; ++j)
//            {
//                dst[i][j] = src[i][j][offset];
//            }
//        }
//    }
//
//
//};

#endif //OPENFPM_PDATA_BLOCKMAPGPU_DIMENSIONALITYWRAPPERS_CUH
