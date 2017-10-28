/*
 * SparseGridUtil.hpp
 *
 *  Created on: Oct 27, 2017
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRIDUTIL_HPP_
#define OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRIDUTIL_HPP_


//! transform T=aggregate<float,double,int> into aggregate<float[n_ele],double[n_ele],int[n_ele]>
template <typename n_ele, typename T>
struct Ft_chunk
{
	typedef std::array<typename std::remove_const<typename std::remove_reference<T>::type>::type,n_ele::value> type;
};

template<unsigned int dim>
struct default_chunking
{
	typedef void type;
};



template<>
struct default_chunking<1>
{
	typedef boost::mpl::vector<boost::mpl::int_<128>> type;

	typedef boost::mpl::vector<boost::mpl::int_<7>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<7>> shift_c;

	typedef boost::mpl::int_<128> size;
};

template<>
struct default_chunking<2>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<8>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<3>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<9>> shift_c;

	typedef boost::mpl::int_<512> size;
};

template<>
struct default_chunking<3>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::int_<1024> size;
};

template<>
struct default_chunking<4>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>,
							   boost::mpl::int_<1>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>,
							   boost::mpl::int_<0>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::int_<1024> size;
};

template<>
struct default_chunking<5>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::int_<1024> size;
};


template<>
struct default_chunking<6>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>> shift;

	typedef boost::mpl::int_<1024> size;
};


template<>
struct default_chunking<7>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::int_<1024> size;
};


template<>
struct default_chunking<8>
{
	typedef boost::mpl::vector<boost::mpl::int_<64>,
			                   boost::mpl::int_<4>,
							   boost::mpl::int_<4>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>,
							   boost::mpl::int_<1>> type;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<2>,
							   boost::mpl::int_<2>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>,
							   boost::mpl::int_<0>> shift;

	typedef boost::mpl::vector<boost::mpl::int_<6>,
			                   boost::mpl::int_<8>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>,
							   boost::mpl::int_<10>> shift_c;

	typedef boost::mpl::int_<1024> size;
};

template<unsigned int dim,
         typename T,
		 typename S,
		 typename layout=typename memory_traits_lin<T>::type,
		 typename chunking = default_chunking<dim>>
class sgrid_cpu
{

};

#define SGRID_CACHE 4


template<unsigned int dim, typename chunk>
struct key_shift
{
	inline static void shift(grid_key_dx<dim> & k)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error dimensionality " << dim << " is not implemented" << std::endl;
	}
};

template<unsigned int dim, typename chunk>
struct sublin
{
	inline static size_t lin(grid_key_dx<dim> & k)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " error dimensionality " << dim << " is not implemented" << std::endl;

		return 0;
	}
};

template<typename chunk>
struct key_shift<1,chunk>
{
	inline static void shift(grid_key_dx<1> & kh, grid_key_dx<1> &kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value -1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
	}

	inline static void cpos(grid_key_dx<1> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
	}
};

template<typename chunk>
struct sublin<1,chunk>
{
	inline static size_t lin(grid_key_dx<1> & k)
	{
		return k.get(0);
	}
};

template<typename chunk>
struct key_shift<2,chunk>
{
	inline static void shift(grid_key_dx<2> & kh, grid_key_dx<2> &kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value-1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value-1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
	}

	inline static void cpos(grid_key_dx<2> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
	}
};

template<typename chunk>
struct sublin<2,chunk>
{
	inline static size_t lin(grid_key_dx<2> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value;
	}
};

template<typename chunk>
struct key_shift<3,chunk>
{
	inline static void shift(grid_key_dx<3> & kh, grid_key_dx<3> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<typename chunk::type,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<typename chunk::shift,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<typename chunk::type,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<typename chunk::shift,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<typename chunk::type,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<typename chunk::shift,boost::mpl::int_<2>>::type::value);
	}

	inline static void cpos(grid_key_dx<3> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<typename chunk::shift,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<typename chunk::shift,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) << boost::mpl::at<typename chunk::shift,boost::mpl::int_<2>>::type::value);
	}
};

template<typename chunk>
struct sublin<3,chunk>
{
	inline static size_t lin(grid_key_dx<3> & k)
	{
		return k.get(0) + (k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value) +
			   (k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
	}
};

template<typename chunk>
struct key_shift<4,chunk>
{
	inline static void shift(grid_key_dx<4> & kh, grid_key_dx<4> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kl.set_d(3,kh.get(3) & (boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value - 1));
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
	}

	inline static void cpos(grid_key_dx<4> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kh.set_d(3,kh.get(3) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
	}
};


template<typename chunk>
struct sublin<4,chunk>
{
	inline static size_t lin(grid_key_dx<4> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value +
			   k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value +
			   k.get(3) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value;
	}
};

template<typename chunk>
struct key_shift<5,chunk>
{
	inline static void shift(grid_key_dx<5> & kh, grid_key_dx<5> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kl.set_d(3,kh.get(3) & (boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value - 1));
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kl.set_d(4,kh.get(4) & (boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value - 1));
		kh.set_d(4,kh.get(4) >> boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
	}

	inline static void cpos(grid_key_dx<5> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kh.set_d(3,kh.get(3) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kh.set_d(4,kh.get(4) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
	}

};

template<typename chunk>
struct sublin<5,chunk>
{
	inline static size_t lin(grid_key_dx<5> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value +
			   k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value +
			   k.get(3) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value +
			   k.get(4) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value;
	}
};

template<typename chunk>
struct key_shift<6,chunk>
{
	inline static void shift(grid_key_dx<6> & kh, grid_key_dx<6> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kl.set_d(3,kh.get(3) & (boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value - 1));
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kl.set_d(4,kh.get(4) & (boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value - 1));
		kh.set_d(4,kh.get(4) >> boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kl.set_d(5,kh.get(5) & (boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value - 1));
		kh.set_d(5,kh.get(5) >> boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
	}

	inline static void cpos(grid_key_dx<6> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kh.set_d(3,kh.get(3) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kh.set_d(4,kh.get(4) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kh.set_d(5,kh.get(5) << boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
	}
};

template<typename chunk>
struct sublin<6,chunk>
{
	inline static size_t lin(grid_key_dx<6> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value +
			   k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value +
			   k.get(3) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value +
			   k.get(4) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value +
			   k.get(5) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value;

	}
};


template<typename chunk>
struct key_shift<7,chunk>
{
	inline static void shift(grid_key_dx<7> & kh, grid_key_dx<7> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kl.set_d(3,kh.get(3) & (boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value - 1));
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kl.set_d(4,kh.get(4) & (boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value - 1));
		kh.set_d(4,kh.get(4) >> boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kl.set_d(5,kh.get(5) & (boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value - 1));
		kh.set_d(5,kh.get(5) >> boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
		kl.set_d(6,kh.get(6) & (boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value - 1));
		kh.set_d(6,kh.get(6) >> boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value);
	}

	inline static void cpos(grid_key_dx<7> & kh)
	{
		kh.set_d(0,kh.get(0) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kh.set_d(3,kh.get(3) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kh.set_d(4,kh.get(4) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kh.set_d(5,kh.get(5) << boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
		kh.set_d(6,kh.get(6) << boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value);
	}
};

template<typename chunk>
struct sublin<7,chunk>
{
	inline static size_t lin(grid_key_dx<7> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value +
			   k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value +
			   k.get(3) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value +
			   k.get(4) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value +
			   k.get(5) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value +
			   k.get(6) << boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value;
	}
};

template<typename chunk>
struct key_shift<8,chunk>
{
	inline static void shift(grid_key_dx<8> & kh, grid_key_dx<8> & kl)
	{
		kl.set_d(0,kh.get(0) & (boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value - 1));
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kl.set_d(1,kh.get(1) & (boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value - 1));
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kl.set_d(2,kh.get(2) & (boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value - 1));
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kl.set_d(3,kh.get(3) & (boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value - 1));
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kl.set_d(4,kh.get(4) & (boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value - 1));
		kh.set_d(4,kh.get(4) >> boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kl.set_d(5,kh.get(5) & (boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value - 1));
		kh.set_d(5,kh.get(5) >> boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
		kl.set_d(6,kh.get(6) & (boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value - 1));
		kh.set_d(6,kh.get(6) >> boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value);
		kl.set_d(7,kh.get(7) & (boost::mpl::at<chunk,boost::mpl::int_<7>>::type::value - 1));
		kh.set_d(7,kh.get(7) >> boost::mpl::at<chunk,boost::mpl::int_<7>>::type::value);
	}

	inline static void cpos(grid_key_dx<8> & kh)
	{
		kh.set_d(0,kh.get(0) >> boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value);
		kh.set_d(1,kh.get(1) >> boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value);
		kh.set_d(2,kh.get(2) >> boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value);
		kh.set_d(3,kh.get(3) >> boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value);
		kh.set_d(4,kh.get(4) >> boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value);
		kh.set_d(5,kh.get(5) >> boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value);
		kh.set_d(6,kh.get(6) >> boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value);
		kh.set_d(7,kh.get(7) >> boost::mpl::at<chunk,boost::mpl::int_<7>>::type::value);
	}
};

template<typename chunk>
struct sublin<8,chunk>
{
	inline static size_t lin(grid_key_dx<8> & k)
	{
		return k.get(0) + k.get(1) << boost::mpl::at<chunk,boost::mpl::int_<0>>::type::value +
			   k.get(2) << boost::mpl::at<chunk,boost::mpl::int_<1>>::type::value +
			   k.get(3) << boost::mpl::at<chunk,boost::mpl::int_<2>>::type::value +
			   k.get(4) << boost::mpl::at<chunk,boost::mpl::int_<3>>::type::value +
			   k.get(5) << boost::mpl::at<chunk,boost::mpl::int_<4>>::type::value +
			   k.get(6) << boost::mpl::at<chunk,boost::mpl::int_<5>>::type::value +
			   k.get(7) << boost::mpl::at<chunk,boost::mpl::int_<6>>::type::value;
	}
};


#endif /* OPENFPM_DATA_SRC_SPARSEGRID_SPARSEGRIDUTIL_HPP_ */
