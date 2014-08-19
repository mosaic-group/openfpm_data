// This class define several reduce type

template <class T>
struct reduce_type
{
	typedef T type;
};

template<> struct reduce_type<float[]>
{
	typedef float type;
};

template<> struct reduce_type<double[]>
{
	typedef double type;
};


template<> struct reduce_type<int[]>
{
	typedef int type;
};
