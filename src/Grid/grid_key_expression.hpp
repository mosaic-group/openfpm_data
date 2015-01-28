#ifndef GRID_KEY_EXPRESSION
#define GRID_KEY_EXPRESSION

/*! \brief Expression template for grid_key_dx
 *
 */

/*! \brief Main class that encapsulate an expression
 *
 * \param dim dimensionality
 * \param exp expression type
 *
 */
template<typename exp>
class grid_key_dx_expression
{
public:

	mem_id value(int i) const
	{
		return static_cast<const exp &>(*this).value(i);
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
class grid_key_dx_sum : public grid_key_dx_expression<grid_key_dx_sum<dim,exp1,exp2>>
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

#endif
