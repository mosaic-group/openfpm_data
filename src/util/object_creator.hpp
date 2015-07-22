#ifndef VECTOR_CREATOR
#define VECTOR_CREATOR

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
 #include <boost/type_traits/remove_reference.hpp>

/*! \brief This is a container for the sending buffers
 *
 * It is used in ghost_get to create a particular object with the properties selected
 *
 * \tparam Is a boost::fusion::vector with the properties selected
 *
 *
 */
template<typename v>
struct object
{
	typedef v type;
	typedef typename boost::fusion::result_of::size<v>::type size_tpy;

	type data;

	static const int max_prop = size_tpy::value;
};

/*! \brief Implementation of vector creator
 *
 * \tparam v boost::fusion::vector
 * \tparam vc result boost::fusion::vector
 *
 */
template<typename v, typename vc, int... prp>
class object_creator_impl
{
};

/*! \brief Implementation of vector creator
 *
 * \tparam v original boost::fusion::vector
 * \tparam vc result boost::fusion::vector
 * \tparam remaining properties to push
 *
 */
template<typename v, typename vc, int p1, int... prp>
struct object_creator_impl<v,vc,p1,prp...>
{
	typedef typename object_creator_impl<v,vc,prp... >::type vc_step;

	typedef typename boost::remove_reference< typename boost::mpl::at< v,boost::mpl::int_<p1> >::type>::type ele;

	// push on the vector the element p1
	typedef typename boost::mpl::push_front<vc_step, ele >::type type;
};

/*! \brief Implementation of vector creator
 *
 * \tparam v original boost::fusion::vector
 * \tparam vc result boost::fusion::vector
 * \tparam remaining properties to push
 */
template<typename v, typename vc, int prp>
struct object_creator_impl<v,vc,prp>
{
	typedef typename boost::remove_reference< typename boost::mpl::at< v,boost::mpl::int_<prp> >::type>::type ele;

	// push on the vector the element p1
	typedef typename boost::mpl::push_front<vc, ele >::type type;
};

/*! \brief It create a boost::fusion vector with the selected properties
 *
 * \tparam v boost::fusion::vector
 * \tparam prp selected properties
 *
 * [Example]
 *
 * typedef boost::fusion::vector<int,float,double,long int, char, short> v;
 *
 * vector_creator<v,0,2,3>::type
 *
 * vector_creator<v,0,2,3>::type is boost::fusion::vector<int,double,long int>
 *
 */

template<typename v, int... prp>
struct object_creator
{
	typedef typename boost::fusion::result_of::as_vector<typename object_creator_impl<v,boost::mpl::vector<>,prp... >::type>::type type;
};

//! specialization when no properties are passed
template<typename v>
struct object_creator<v>
{
	typedef typename boost::fusion::vector<> type;
};

#endif
