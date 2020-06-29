/*
 * map_vector_printers.hpp
 *
 *  Created on: Feb 10, 2019
 *      Author: i-bird
 */

#ifndef MAP_VECTOR_PRINTERS_HPP_
#define MAP_VECTOR_PRINTERS_HPP_


/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to print the elements of the vector
 *
 * \tparam encap source
 * \tparam encap dst
 *
 */
template<typename vector_type, unsigned int ... prp>
struct vector_printer
{
	//! element to print
	size_t & ele;

	//! vector to print
	vector_type & vt;

	// stringstream
	std::stringstream & ss;

	typedef typename to_boost_vmpl<prp...>::type vprp;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline vector_printer(vector_type & vt, size_t & ele, std::stringstream & ss)
	:vt(vt),ele(ele),ss(ss)
	{};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		ss << vt.template get<T::value>(ele) << "    ";
	}
};


#endif /* MAP_VECTOR_PRINTERS_HPP_ */
