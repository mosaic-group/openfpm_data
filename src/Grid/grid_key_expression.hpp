#ifndef GRID_KEY_EXPRESSION
#define GRID_KEY_EXPRESSION

#include "util/common.hpp"

template<unsigned int dim, typename index_type = long int> class grid_key_dx;
template<int dim, typename exp1, typename exp2> class grid_key_dx_sum;
template<int dim, typename exp1, typename exp2> class grid_key_dx_sub;


/*! \brief Expression template for grid_key_dx
 *
 */

/*! \brief Main class that encapsulate an expression
 *
 * \param dim dimensionality
 * \param exp expression type
 *
 */
template<unsigned int dim, typename exp>
class grid_key_dx_expression
{
public:

	mem_id value(int i) const
	{
		return static_cast<const exp &>(*this).value(i);
	}

	/* \brief subtract this expression with another expression
	 *
	 * \param key to subtract
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	inline grid_key_dx_sub<dim,grid_key_dx_expression<dim,exp>,grid_key_dx<dim>> operator-(const grid_key_dx<dim> & key) const
	{
		grid_key_dx_sub<dim,grid_key_dx_expression<dim,exp>,grid_key_dx<dim>> exp_sum(*this,key);

		return exp_sum;
	}

	/* \brief subtract this expression a grid key
	 *
	 * \param key to subtract
	 *
	 * \return a grid_key_dx_expression that encapsulate the expression
	 *
	 */
	template <typename T> inline grid_key_dx_sub<dim,grid_key_dx_expression<dim,exp>,grid_key_dx_expression<dim,T> > operator-(const grid_key_dx_expression<dim,T> & key) const
	{
		grid_key_dx_sub< dim,grid_key_dx_expression<dim,exp>,grid_key_dx_expression<dim,T> > exp_sum(*this,key);

		return exp_sum;
	}
};


/*! \brief Main class that encapsulate a sum expression
 *
 * \param dim dimensionality
 * \param exp1 expression 1
 * \param exp2 expression 2
 *
 */
template<int dim, typename exp1, typename exp2>
class grid_key_dx_sum : public grid_key_dx_expression<dim,grid_key_dx_sum<dim,exp1,exp2>>
{
	const exp1 & e1;
	const exp2 & e2;

public:

	grid_key_dx_sum(const exp1 & ex1, const exp2 & ex2)
	:e1(ex1),e2(ex2)
	{}

	mem_id value(int i) const
	{
		return e1.value(i) + e2.value(i);
	}
};

/*! \brief Main class that encapsulate a sub expression
 *
 * \param dim dimensionality
 * \param exp1 expression 1
 * \param exp2 expression 2
 *
 */
template<int dim, typename exp1, typename exp2>
class grid_key_dx_sub : public grid_key_dx_expression<dim,grid_key_dx_sub<dim,exp1,exp2>>
{
	const exp1 & e1;
	const exp2 & e2;

public:

	grid_key_dx_sub(const exp1 & ex1, const exp2 & ex2)
	:e1(ex1),e2(ex2)
	{}

	mem_id value(int i) const
	{
		return e1.value(i) - e2.value(i);
	}
};

#endif
