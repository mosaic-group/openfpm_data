// This class define several reduce type

//! In general a reduction of a type T produce a type T
template <class T>
struct reduce_type
{
	typedef T type;
};

//! A reduction operation on an array of float is a float
template<> struct reduce_type<float[]>
{
	//! define the float
	typedef float type;
};

//! A reduction operation on an array of double is a double
template<> struct reduce_type<double[]>
{
	//! define the double
	typedef double type;
};

//! A reduction operation on an array of int is an int
template<> struct reduce_type<int[]>
{
	//! define the int
	typedef int type;
};
