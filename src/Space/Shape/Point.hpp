#ifndef POINT_HPP
#define POINT_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "Grid/Encap.hpp"
#include "Point_operators.hpp"


/*! \brief This class implement the point shape in an N-dimensional space
 *
 * \param T type of the space
 * \param dim dimensionality
 *
 */

template<unsigned int dim ,typename T> class Point
{
	public:

	//! Indicate that this object is vtk writable
	typedef int is_vtk_writable;

	//! coordinate type
	typedef T coord_type;

	//! boost fusion that store the point
	typedef boost::fusion::vector<T[dim]> type;

	//! Indicate that is a Point
	typedef int yes_is_point;

	//! structure that store the data of the point
	type data;

	//! Property id of the point
	static const unsigned int x = 0;


	/*! \brief Evaluate the expression and save the result on the point
	 *
	 * \param p_exp point expression to evaluate
	 *
	 */
	template<typename orig, typename exp1, typename exp2, unsigned int op> __device__ __host__ Point(const point_expression_op<orig,exp1,exp2,op> & p_exp)
	{
		this->operator=(p_exp);
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	__device__ __host__ inline Point(const Point<dim,T> && p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Point constructor from point
	 *
	 * \param p the point
	 *
	 */
	__device__ __host__ inline Point(const Point<dim,T> & p)
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p.get(i);}
	}

	/*! \brief Constructor from an array
	 *
	 * \param p array with the coordinate of the point
	 *
	 */
	__device__ __host__ inline Point(const T (&p)[dim])
	{
	    for(size_t i = 0; i < dim ; i++)
	    {get(i) = p[i];}
	}

	/*! \brief Constructor from scalar
	 *
	 * \param d scalar
	 *
	 */
	__device__ __host__ inline Point(T d)
	{
		this->operator=(d);
	}

	/*! \brief Point constructor
	 *
	 * \param p Point
	 *
	 */
	template <typename S> inline Point(const Point<dim,S> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{get(i) = static_cast<S>(p.get(i));}
	}

	/*! \brief Point constructor
	 *
	 * \param p encapc Point
	 *
	 */
	template <unsigned int d, typename M> inline Point(const encapc<d,Point<dim,T>,M> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) = p.template get<0>()[i];
	}

	/*! \brief Point constructor from multi array
	 *
	 *
	 *
	 */
	template <typename vmpl> inline __device__ __host__  Point(const openfpm::detail::multi_array::sub_array_openfpm<T,1,vmpl> & mar)
	{
		for (unsigned int i = 0 ; i < dim ; i++)
		{get(i) = mar[i];}
	}

	/*! \brief Point constructor from multi array
	 *
	 *
	 *
	 */
	template <typename vmpl> inline __device__ __host__  Point(const openfpm::detail::multi_array::const_sub_array_openfpm<T,1,vmpl> & mar)
	{
		for (unsigned int i = 0 ; i < dim ; i++)
		{get(i) = mar[i];}
	}

	/*! \brief Constructor from a list
	 *
	 * [Example] Point<3,float> p({0.0,0.0,1.0})
	 *
	 * \param p1 initializer list
	 *
	 */
	__device__ __host__ inline Point(std::initializer_list<T> p1)
	{
		size_t i = 0;
	    for(T x : p1)
	    {get(i) = x;i++;}
	}

	//! Default contructor
	__device__ __host__  inline Point()
	{}

	/*! \brief Get coordinate
	 *
	 * \param i dimension
	 * \return the i-coordinate of the point
	 *
	 */
	__device__ __host__ inline const T & get(int i) const
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Get coordinate
	 *
	 * \param i dimension
	 * \return the i-coordinate of the point
	 *
	 */

	inline T get_vtk(size_t i) const
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Get coordinate
	 *
	 * \param i dimension
	 * \return the i-coordinate of the point
	 *
	 */
	__device__ __host__ inline T& get(int i)
	{
		return boost::fusion::at_c<x>(data)[i];
	}

	/*! \brief Get the component i
	 *
	 * \param i component
	 *
	 * \return the i-component
	 *
	 */

	inline T& operator[](size_t i)
	{
		return get(i);
	}

	/*! \brief Get the component i
	 *
	 * \param i component
	 *
	 * \return the i-component
	 *
	 */

	inline const T& operator[](size_t i) const
	{
		return get(i);
	}

	/*! \brief norm of the vector
	 *
	 * \return the norm of the vector
	 *
	 */
	__device__ __host__ T norm()
	{
		T n = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			n+=get(i) * get(i);

		return sqrt(n);
	}

	/*! \brief  It calculate the distance between 2 points
	 *
	 * The distance between itself (p) and the other point (q)
	 *
	 * \param q target point
	 *
	 * \return the distance
	 *
	 */
	__device__ __host__ T distance(const Point<dim,T> & q) const
	{
		T tot = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			tot += (this->get(i)  - q.get(i)) * (this->get(i)  - q.get(i));

		return sqrt(tot);
	}

	/*! \brief  It calculate the square distance between 2 points
	 *
	 * The distance between itself (p) and the other point (q)
	 *
	 * \param q target point
	 *
	 * \return the square of the distance
	 *
	 */
	T distance2(const Point<dim,T> & q) const
	{
		T tot = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			tot += (this->get(i)  - q.get(i)) * (this->get(i)  - q.get(i));

		return tot;
	}


	/*! \brief Set to zero the point coordinate
	 *
	 *
	 */
	inline void zero()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) = 0;
		}
	}

	/*! \brief Set to one the point coordinate
	 *
	 *
	 */
	inline void one()
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			get(i) = 1;
		}
	}

	/*! \brief Create a point set to zero
	 *
	 * \return a point with all coorfinate set to 0
	 *
	 */
	inline static Point<dim,T> zero_p()
	{
		Point<dim,T> p;

		for (size_t i = 0 ; i < dim ; i++)
		{
			p.get(i) = 0;
		}

		return p;
	}

	/*! \brief Convert the point into a string
	 *
	 * \return the string
	 *
	 */
	std::string toPointString() const
	{
		std::stringstream ps;

		for (size_t i = 0 ; i < dim ; i++)
			ps << "x[" << i << "]=" << get(i) << " ";

		ps << "\n";

		return ps.str();
	}

	/*! \brief exchange the data of two points
	 *
	 * \param p Point to swap with
	 *
	 */
	void swap(Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			T tmp = get(i);
			get(i) = p.get(i);
			p.get(i) = tmp;
		}
	}

	/*! \brief Check if two points match
	 *
	 * \param p point to compare with
	 *
	 * \return true if two points match
	 *
	 */
	__device__ __host__ inline bool operator==(const Point<dim,T> & p) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (p.get(i) != get(i))
				return false;
		}

		return true;
	}

	/*! \brief Check if two points match
	 *
	 * \param p point to compare with
	 *
	 * \return true if two points match
	 *
	 */
	__device__ __host__ inline bool operator!=(const Point<dim,T> & p) const
	{
		return !this->operator==(p);
	}

	/*! \brief Return the string with the point coordinate
	 *
	 * \return the string
	 *
	 */
	std::string to_string() const
	{
		return toString();
	}

	/*! \brief Return the string with the point coordinate
	 *
	 * \return the string
	 *
	 */
	std::string toString() const
	{
		std::string str;

		for (size_t i = 0 ; i < dim - 1 ; i++)
		{
			/* coverty[dead_error_line] */
			str += std::to_string(static_cast<double>(get(i))) + " ";
		}
		str += std::to_string(static_cast<double>(get(dim-1)));

		return str;
	}

	/*! \brief Return the reference to the value at coordinate i
	 *
	 * \param i coordinate to return
	 *
	 * \return the reference
	 *
	 */
	T & value(size_t i)
	{
		return get(i);
	}

	/*! \brief Return the value at coordinate i
	 *
	 * \param i coordinate to return
	 *
	 * \return the value
	 *
	 */
	inline T value(size_t i) const
	{
		return get(i);
	}

	/*! \brief Return the coordinated of the point as reference array
	 *
	 * \return the reference array
	 *
	 */
	T (&asArray())[dim]
	{
		return boost::fusion::at_c<x>(data);
	}

	/*! Convert the point from Point<dim,T> to Point<dim,A>
	 *
	 * \return the converted point
	 *
	 */
	template<typename A> Point<dim,A> convertPoint() const
	{
		Point<dim,A> p;

		for (size_t i = 0; i < dim ; i++)
			p.get(i) = static_cast<A>(get(i));

		return p;
	}


	//! This structure has no internal pointers
	static bool noPointers()
	{
		return true;
	}

	////////////////////////////////////////////////////////////////
	////////////////////// ARITMETIC OPERATORS /////////////////////
	////////////////////////////////////////////////////////////////

	/*! \brief Fill the vector with the evaluated expression
	 *
	 * \param p_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<typename orig, typename exp1, typename exp2, unsigned int op> __device__ __host__ Point<dim,T> & operator=(const point_expression_op<orig,exp1,exp2,op> & p_exp)
	{
		p_exp.init();

		for (size_t i = 0; i < dim ; i++)
			get(i) = p_exp.value(i);

		return *this;
	}

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \param p_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ Point<dim,T> & operator=(const point_expression<T[dim]> & p_exp)
	{
		p_exp.init();

		for (size_t i = 0; i < dim ; i++)
			get(i) = p_exp.value(i);

		return *this;
	}

	/*! \brief Fill the point with the value specified in the array
	 *
	 * \param p array
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ Point<dim,T> & operator=(const T (& p)[dim])
	{
		for (size_t i = 0; i < dim ; i++)
			get(i) = p[i];

		return *this;
	}

	/*! \brief Fill the vector property with the evaluated expression
	 *
	 * \tparam check disable this method if T is a constant
	 *
	 * \param p_exp expression to evaluate
	 *
	 * \return itself
	 *
	 */
	template<typename T1, typename check = typename std::enable_if<std::is_const<T1>::value == false>::type> __device__ __host__ Point<dim,T> & operator=(const point_expression<const T1[dim]> & p_exp)
	{
		p_exp.init();

		for (size_t i = 0; i < dim ; i++)
			get(i) = p_exp.value(i);

		return *this;
	}

	/*! \brief divide each component by an array
	 *
	 * \param ar Component wise division
	 *
	 * \return itself
	 *
	 */
	template<typename aT> __device__ __host__ inline Point<dim,T> operator/(const aT (&ar)[dim])
	{
		Point<dim,T> result;

		for (size_t i = 0 ; i < dim ; i++)
			result.get(i) = get(i) / ar[i];

		return result;
	}

	/*! \brief divide each component by a constant
	 *
	 * \param c constant
	 *
	 * \return itself
	 *
	 */
	template<typename aT> __device__ __host__ inline Point<dim,T> operator/=(const aT c)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) = get(i) / c;

		return *this;
	}


	/*! \brief Fill the vector property with some value
	 *
	 * \param d value to fill
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ Point<dim,T> & operator=(T d)
	{
		for (size_t i = 0; i < dim ; i++)
			get(i) = d;

		return *this;
	}

	/*! \brief operator= between points
	 *
	 * \param p Point
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline Point<dim,T> & operator=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) = p.get(i);

		return *this;
	}

	/*! \brief Subtract two points
	 *
	 * \param p point to subtract
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline Point<dim,T> & operator-=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
			get(i) -= p.get(i);

		return *this;
	}

	/*! \brief Sum two points
	 *
	 * \param p point to sum
	 *
	 * \return itself
	 *
	 */
	__device__ __host__ inline Point<dim,T> & operator+=(const Point<dim,T> & p)
	{
		for (size_t i = 0 ; i < dim ; i++)
		{get(i) += p.get(i);}

		return *this;
	}

	/*! \brief Do nothing stub operation
	 *
	 * Required to make the code compilable
	 *
	 */
	__device__ __host__ inline void init() const
	{}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	//! The point has one property
	static const unsigned int max_prop = 1;
	static const unsigned int max_prop_real = 1;

	//! expose the dimension
	static const unsigned int dims = dim;

	//! expose the dimension with a different name
	static const unsigned int nvals = dim;
};

/*! \brief Convert an array of point coordinate into string
 *
 * \param p coordinate on each dimension
 *
 * \return the string
 *
 */
template <unsigned int N, typename T> std::string toPointString(const T (&p)[N] )
{
	std::stringstream ps;

	for (size_t i = 0 ; i < N ; i++)
		ps << "x[" << i << "]=" << p[i] << " ";

	ps << "\n";

	return ps.str();
}

/*! \brief Convert an encapsulated point into string
 *
 * \param p coordinate on each dimension
 *
 * \return the string
 *
 */
template <unsigned int N, typename T, typename Mem> std::string toPointString(const encapc<1,Point<N,T>,Mem> & p )
{
	std::stringstream ps;

	for (size_t i = 0 ; i < N ; i++)
		ps << "x[" << i << "]=" << p.template get<Point<N,T>::x>()[i] << " ";

	ps << "\n";

	return ps.str();
}


//! A point is a vector on a computer (But do not say this to a Mathematician)

template<unsigned int dim, typename T>  using VectorS = Point<dim,T>;

template<typename T, typename Sfinae = void>
struct is_Point: std::false_type {};


/*! \brief Check if a type T is an aggregate
 *
 * return true if T is an aggregate
 *
 */
template<typename T>
struct is_Point<T, typename Void< typename T::yes_is_point>::type> : std::true_type
{};

#endif
