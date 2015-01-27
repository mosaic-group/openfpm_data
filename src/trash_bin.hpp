/*
 * trash_bin.hpp
 *
 *  Created on: Nov 24, 2014
 *      Author: i-bird
 */

#ifndef TRASH_BIN_HPP_
#define TRASH_BIN_HPP_


/*! \brief Edge class that encapsulate an object T
 *
 * Edge class that encapsulate an object T that store the property of the edge
 * It is basically used to append a field size_t
 * to the original object T, this field is the neighborhood vertex in the adjacency list
 *
 */

template<typename T>
class E_p
{
public:
	//! insert a size_t property to the boost::vector
	typedef typename boost::fusion::result_of::push_front<typename T::type, size_t>::type E_a;

	//! type for the edge
	typedef E_a type;

	E_a data;
};

/*! \brief Vertex class that encapsulate an object T
 *
 * Vertex class that encapsulate an object T that store the property of the vertex
 * It is basically used to insert a field size_t
 * to the original object T, this field is the end list point in the adjacency list
 *
 */

template<typename T>
class V_p
{
public:
	//! insert a size_t property to the boost::vector
	typedef typename boost::fusion::result_of::push_front<typename T::type, size_t>::type V_a;

	//! type for the vertex
	typedef V_a type;

	V_a data;
};


#endif /* TRASH_BIN_HPP_ */


/*! \brief Get the reference of the selected element
 *
 * Get the reference of the selected element
 *
 * \param p property to get (is an integer)
 * \param v1 grid_key that identify the element in the grid
 *
 */
/*		template <unsigned int p>inline typename type_cpu_prop<p,T>::type & getBoostVector(grid_key_dx<dim> & v1)
{
#ifdef MEMLEAK_CHECK
	check_valid(&boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1))));
#endif
	return boost::fusion::at_c<p>(data.mem_r->operator[](g1.LinId(v1)));
}*/
