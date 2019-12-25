#ifndef POINT_TEST_HPP
#define POINT_TEST_HPP

#include <boost/fusion/sequence/intrinsic/at_c.hpp>
#include <boost/fusion/include/at_c.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>
#include "boost/multi_array.hpp"
#include "Point_orig.hpp"
#include "Grid/Encap.hpp"
#include "data_type/aggregate.hpp"


/*! \brief Test structure used for several test
 *
 * It is a test structure used for several test it contain 4 scalar "x,y,z,s"
 * one vector property v[3] and one tensor or rank 2 t[3][3]
 *
 * It is the OpenFPM data structure format for type parsing of see openFPM_data wiki
 * for more information
 *
 * ### Declaration of a point
 * \snippet vector_test_util.hpp Point declaration
 *
 * ### Create a type definition
 *
 * \snippet vector_test_util.hpp typedef point
 *
 * ### Access the members
 * \snippet vector_test_util.hpp Point usage
 *
 */
template<typename T> class Point2D_test
{
public:

	//! declaration of what the Point2D_test store
	typedef boost::fusion::vector<T,T,T,T,T[2],T[2][2]> type;

	//! declaration of what the Point2D_test store
	typedef boost::fusion::vector<T,T,T,T,T[2],T[2][2]> type_real;

	//! in case of usage with staggered grid it define which properties are staggered in the cell grid
	static constexpr bool stag_mask[]={false,false,false,true,true,true};

	//! object itself
	type data;

	//! x property is at position 0 in the boost::fusion::vector
	static const unsigned int x = 0;

	//! y property is at position 1 in the boost::fusion::vector
	static const unsigned int y = 1;

	//! z property is at position 2 in the boost::fusion::vector
	static const unsigned int z = 2;

	//! s property is at position 3 in the boost::fusion::vector
	static const unsigned int s = 3;

	//! v property is at position 4 in the boost::fusion::vector
	static const unsigned int v = 4;

	//! t property is at position 5 in the boost::fusion::vector
	static const unsigned int t = 5;

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop = 6;

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop_real = 6;

	/*! \brief indicate that this structure has no pointers
	 *
	 * \return true
	 *
	 */
	static inline bool noPointers()
	{
		return true;
	}
};

/*! \brief Test structure used for several test
 *
 * It is a test structure used for several test it contain 4 scalar "x,y,z,s"
 * one vector property v[3] and one tensor or rank 2 t[3][3]
 *
 * It is the format for type parsing of in the openfpm structures see openFPM_data wiki
 * for more information
 *
 * ### Declaration of a point
 * \snippet vector_test_util.hpp Point declaration
 *
 * ### Create a type definition
 *
 * \snippet vector_test_util.hpp typedef point
 *
 * ### Access the members
 * \snippet vector_test_util.hpp Point usage
 *
 */
template<typename T> class Point_test
{
public:

#ifdef SE_CLASS3

	//! declaration of what the Point_test store
	typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3],SE3_ADD_PROP(6)> type;

	//! declaration of what the Point_test store
	typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type_real;

#else

	//! declaration of what the Point_test store
	typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type;

	//! declaration of what the Point_test store
	typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type_real;

#endif

	//! in case usage with a staggered grid indicate which properties are staggered in the cell
	static constexpr bool stag_mask[]={false,false,false,true,true,true};

	//! The object itself
	type data;

	//! x property is at position 0 in the boost::fusion::vector
	static const unsigned int x = 0;

	//! y property is at position 1 in the boost::fusion::vector
	static const unsigned int y = 1;

	//! z property is at position 2 in the boost::fusion::vector
	static const unsigned int z = 2;

	//! s property is at position 3 in the boost::fusion::vector
	static const unsigned int s = 3;

	//! v property is at position 4 in the boost::fusion::vector
	static const unsigned int v = 4;

	//! t property is at position 5 in the boost::fusion::vector
	static const unsigned int t = 5;

#ifdef SE_CLASS3

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop = SE3_MAX_PROP(6);

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop_real = 6;

#else

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop = 6;

	//! number of properties in the boost::fusion::vector
	static const unsigned int max_prop_real = 6;

#endif

	// Setter method

	/*! \brief set the x property
	 *
	 * \param x_
	 *
	 */
	inline void setx(T x_)	{boost::fusion::at_c<0>(data) = x_;};

	/*! \brief set the y property
	 *
	 * \param y_
	 *
	 */
	inline void sety(T y_)	{boost::fusion::at_c<1>(data) = y_;};

	/*! \brief set the z property
	 *
	 * \param z_
	 *
	 */
	inline void setz(T z_)	{boost::fusion::at_c<2>(data) = z_;};

	/*! \brief set the s property
	 *
	 * \param s_
	 *
	 */
	inline void sets(T s_)	{boost::fusion::at_c<3>(data) = s_;};

	/*! \brief set the v property
	 *
	 * \param i component to set
	 * \param v_ value
	 *
	 */
	inline void setv(size_t i,T v_)	{boost::fusion::at_c<4>(data)[i] = v_;}

	/*! \brief set the t property
	 *
	 * \param i component to set
	 * \param j component to set
	 * \param t_ value
	 *
	 */
	inline void sett(size_t i, size_t j,T t_)	{boost::fusion::at_c<5>(data)[i][j] = t_;}


	//! getter method for a general property i
	template<unsigned int i>
	inline auto get() -> decltype(boost::fusion::at_c<i>(data))
	{return boost::fusion::at_c<i>(data);}

	//! getter method for a general property i
	template<unsigned int i>
	inline auto get() const -> decltype(boost::fusion::at_c<i>(data))
	{return boost::fusion::at_c<i>(data);}

	//! Default constructor
	Point_test()
	{}

	/*! \brief check if two point match
	 *
	 * \param p point to compare
	 *
	 */
	bool operator==(const Point_test<float> & p) const
	{
		if (boost::fusion::at_c<0>(data) != boost::fusion::at_c<0>(p.data))	return false;
		if (boost::fusion::at_c<1>(data) != boost::fusion::at_c<1>(p.data))	return false;
		if (boost::fusion::at_c<2>(data) != boost::fusion::at_c<2>(p.data))	return false;
		if (boost::fusion::at_c<3>(data) != boost::fusion::at_c<3>(p.data))	return false;

		for (size_t i = 0 ; i < 3 ; i++)
			if (boost::fusion::at_c<4>(data)[i] != boost::fusion::at_c<4>(p.data)[i])	return false;

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				if (boost::fusion::at_c<5>(data)[i][j] != boost::fusion::at_c<5>(p.data)[i][j])	return false;
			}
		}

		return true;
	}

	/*! \brief Sum the point
	 *
	 * \param p point to sum
	 *
	 * \return this
	 *
	 */
	Point_test<float> & operator+=(const Point_test<float> & p)
	{
		boost::fusion::at_c<0>(data) += boost::fusion::at_c<0>(p.data);
		boost::fusion::at_c<1>(data) += boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) += boost::fusion::at_c<2>(p.data);
		boost::fusion::at_c<3>(data) += boost::fusion::at_c<3>(p.data);

		for (size_t i = 0 ; i < 3 ; i++)
			boost::fusion::at_c<4>(data)[i] += boost::fusion::at_c<4>(p.data)[i];

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
				boost::fusion::at_c<5>(data)[i][j] += boost::fusion::at_c<5>(p.data)[i][j];
		}

		return *this;
	}

	/*! \brief Copy constructor from encapc (encapsulated point)
	 *
	 * \param p ecapsulated point
	 *
	 */
	template <unsigned int dim, typename Mem> inline Point_test(const encapc<dim,Point_test<T>,Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();
		boost::fusion::at_c<3>(data) = p.template get<3>();

		for (size_t i = 0 ; i < 3 ; i++)
			boost::fusion::at_c<4>(data)[i] = p.template get<4>()[i];

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				boost::fusion::at_c<5>(data)[i][j] = p.template get<5>()[i][j];
			}
		}
	}

	/*! \brief constructor from another point
	 *
	 *  \param p point to copy
	 *
	 */
	inline Point_test(const Point_test<T> & p)
	{
		boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
		boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
		boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);

		for (size_t i = 0 ; i < 3 ; i++)
			boost::fusion::at_c<4>(data)[i] = boost::fusion::at_c<4>(p.data)[i];

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				boost::fusion::at_c<5>(data)[i][j] = boost::fusion::at_c<5>(p.data)[i][j];
			}
		}
	}

	/*! \brief Copy the point
	 *
	 * \param p point
	 *
	 * \return this
	 *
	 */
	__host__ __device__ inline Point_test<T> operator= (const Point_test<T> & p)
	{
		boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
		boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
		boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);

		for (size_t i = 0 ; i < 3 ; i++)
			boost::fusion::at_c<4>(data)[i] = boost::fusion::at_c<4>(p.data)[i];

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				boost::fusion::at_c<5>(data)[i][j] = boost::fusion::at_c<5>(p.data)[i][j];
			}
		}

		return *this;
	}

	/*! \brief noPointers function
	 *
	 * It notify that Point_test does not have any pointer and is safe to send
	 *
	 * \return true
	 *
	 */
	static bool noPointers()	{return true;}

	/*! \brief fill
	 *
	 * Fill the point with data
	 *
	 */
	void fill()
	{
		boost::fusion::at_c<0>(data) = 1;
		boost::fusion::at_c<1>(data) = 2;
		boost::fusion::at_c<2>(data) = 3;
		boost::fusion::at_c<3>(data) = 4;

		for (size_t i = 0 ; i < 3 ; i++)
			boost::fusion::at_c<4>(data)[i] = 5;

		for (size_t i = 0 ; i < 3 ; i++)
		{
			for (size_t j = 0 ; j < 3 ; j++)
			{
				boost::fusion::at_c<5>(data)[i][j] = 6;
			}
		}
	}
};


/*! \brief Test structure used for several test
 *
 * It is a test structure used for several test it contain 4 scalar "x,y,z,s"
 * one vector property v[3] and one tensor or rank 2 t[3][3] + the definition
 * of properties names
 *
 * It is the format for type parsing of in the openfpm structures see openFPM_data wiki
 * for more information
 *
 * ### Declaration of a point
 * \snippet vector_test_util.hpp Point prp declaration
 *
 * ### Create a type definition
 *
 * \snippet vector_test_util.hpp typedef point
 *
 * ### Access the members
 * \snippet vector_test_util.hpp Point prp usage
 *
 */
template<typename T> class Point_test_prp
{
public:

  //! declaration of what the Point_test_prp store
  typedef boost::fusion::vector<T,T,T,T,T[3],T[3][3]> type;

  //! Object itself
  type data;

  //! x property is at position 0 in the boost::fusion::vector
  static const unsigned int x = 0;

  //! y property is at position 1 in the boost::fusion::vector
  static const unsigned int y = 1;

  //! z property is at position 2 in the boost::fusion::vector
  static const unsigned int z = 2;

  //! s property is at position 3 in the boost::fusion::vector
  static const unsigned int s = 3;

  //! v property is at position 4 in the boost::fusion::vector
  static const unsigned int v = 4;

  //! t property is at position 5 in the boost::fusion::vector
  static const unsigned int t = 5;

  //! maximum number of properties
  static const unsigned int max_prop = 6;

  // Setter method

  /*! \brief set the x property
   *
   * \param x_
   *
   */
  inline void setx(T x_)	{boost::fusion::at_c<0>(data) = x_;};

  /*! \brief set the y property
   *
   * \param y_
   *
   */
  inline void sety(T y_)	{boost::fusion::at_c<1>(data) = y_;};

  /*! \brief set the z property
   *
   * \param z_
   *
   */
  inline void setz(T z_)	{boost::fusion::at_c<2>(data) = z_;};

  /*! \brief set the s property
   *
   * \param s_
   *
   */
  inline void sets(T s_)	{boost::fusion::at_c<3>(data) = s_;};

  //! Attributes name
  struct attributes
  {
	//! array of attributes name
    static const std::string name[];
  };

  //! getter method for a general property i
  template<unsigned int i> inline typename boost::fusion::result_of::at<type, boost::mpl::int_<i> >::type get()	{return boost::fusion::at_c<i>(data);}

  //! Default constructor
  Point_test_prp()
  {}

  //! constructor from encapc
  template <typename Mem> inline Point_test_prp(const encapc<1,Point_test_prp<T>,Mem> & p)
  {
	  boost::fusion::at_c<0>(data) = p.template get<0>();
	  boost::fusion::at_c<1>(data) = p.template get<1>();
	  boost::fusion::at_c<2>(data) = p.template get<2>();
	  boost::fusion::at_c<3>(data) = p.template get<3>();

	  for (size_t i = 0 ; i < 3 ; i++)
		  boost::fusion::at_c<4>(data)[i] = p.template get<4>()[i];

	  for (size_t i = 0 ; i < 3 ; i++)
	  {
		  for (size_t j = 0 ; j < 3 ; j++)
		  {
			  boost::fusion::at_c<5>(data)[i][j] = p.template get<5>()[i][j];
		  }
	  }
  }

  //! constructor from another point
  inline Point_test_prp(const Point_test_prp<T> & p)
  {
	  boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
	  boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
	  boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
	  boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);

	  for (size_t i = 0 ; i < 3 ; i++)
		  boost::fusion::at_c<4>(data)[i] = boost::fusion::at_c<4>(p.data)[i];

	  for (size_t i = 0 ; i < 3 ; i++)
	  {
		  for (size_t j = 0 ; j < 3 ; j++)
		  {
			  boost::fusion::at_c<5>(data)[i][j] = boost::fusion::at_c<5>(p.data)[i][j];
		  }
	  }
  }


	/*! \brief Copy the point
	 *
	 * \param p point
	 *
	 * \return this
	 *
	 */
  inline Point_test_prp<T> operator= (const Point_test<T> & p)
  {
	  boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
	  boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
	  boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
	  boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);

	  for (size_t i = 0 ; i < 3 ; i++)
		  boost::fusion::at_c<4>(data)[i] = boost::fusion::at_c<4>(p.data)[i];

	  for (size_t i = 0 ; i < 3 ; i++)
	  {
		  for (size_t j = 0 ; j < 3 ; j++)
		  {
			  boost::fusion::at_c<5>(data)[i][j] = boost::fusion::at_c<5>(p.data)[i][j];
		  }
	  }

	  return *this;
  }

	static inline bool noPointers()
	{
		return true;
	}
};

template<typename T> const std::string Point_test_prp<T>::attributes::name[] = {"x","y","z","s","v","t"};

//! point test with only scalar properties
template<typename T> class Point_test_scal
{
public:

  //! declaration of what the Point_test_scal store
  typedef boost::fusion::vector<T,T,T,T> type;

  //! The data itself
  type data;

  //! x property is at position 0 in the boost::fusion::vector
  static const unsigned int x = 0;

  //! y property is at position 1 in the boost::fusion::vector
  static const unsigned int y = 1;

  //! z property is at position 0 in the boost::fusion::vector
  static const unsigned int z = 2;

  //! s property is at position 0 in the boost::fusion::vector
  static const unsigned int s = 3;

  //! the number of properties
  static const unsigned int max_prop = 4;

  // Setter method

  //! set the property x
  inline void setx(T x_)	{boost::fusion::at_c<0>(data) = x_;};
  //! set the property y
  inline void sety(T y_)	{boost::fusion::at_c<1>(data) = y_;};
  //! set the property z
  inline void setz(T z_)	{boost::fusion::at_c<2>(data) = z_;};
  //! set the property s
  inline void sets(T s_)	{boost::fusion::at_c<3>(data) = s_;};

  //! Attributes name
  struct attributes
  {
    //! array of names
    static const std::string name[];
  };

  //! getter method for the property i
  template<unsigned int i> inline typename boost::fusion::result_of::at<type, boost::mpl::int_<i> >::type get()	{return boost::fusion::at_c<i>(data);}

  //! Default constructor
  Point_test_scal()
  {}

  //! constructor from encapc
  template <typename Mem> inline Point_test_scal(const encapc<1,Point_test_scal<T>,Mem> & p)
  {
	  boost::fusion::at_c<0>(data) = p.template get<0>();
	  boost::fusion::at_c<1>(data) = p.template get<1>();
	  boost::fusion::at_c<2>(data) = p.template get<2>();
	  boost::fusion::at_c<3>(data) = p.template get<3>();
  }

  //! constructor from another point
  inline Point_test_scal(const Point_test_scal<T> & p)
  {
	  boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
	  boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
	  boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
	  boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);
  }

  //! operator=
  inline Point_test_scal<T> operator= (const Point_test_scal<T> & p)
  {
	  boost::fusion::at_c<0>(data) = boost::fusion::at_c<0>(p.data);
	  boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
	  boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
	  boost::fusion::at_c<3>(data) = boost::fusion::at_c<3>(p.data);

	  return *this;
  }

	static inline bool noPointers()
	{
		return true;
	}
};

template<typename T> const std::string Point_test_scal<T>::attributes::name[] = {"x","y","z","s"};

#endif
