/*
 * Point_operators_functions.hpp
 *
 *  Created on: Jun 20, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_FUNCTIONS_HPP_
#define OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_FUNCTIONS_HPP_

#define CREATE_ARG_FUNC_CLASS(fun_base,fun_name,OP_ID) \
\
template <typename orig, typename exp1, typename exp2>\
class point_expression_op<orig,exp1,exp2, OP_ID >\
{\
	const exp1 o1;\
\
	mutable typename std::remove_const<typename orig::coord_type>::type scal;\
\
public:\
\
	typedef orig orig_type;\
\
	typedef int is_expression;\
\
	typedef int has_init;\
\
	typedef typename orig::coord_type return_type;\
\
	static const unsigned int nvals = exp1::nvals;\
\
	inline point_expression_op(const exp1 & o1)\
	:o1(o1),scal(0)\
	{}\
\
	inline void init() const\
	{\
		o1.init();\
	}\
\
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > __device__ __host__  inline r_type value(size_t k) const\
	{\
		return fun_base(o1.value(k));\
	}\
\
	template <typename T> __device__ __host__  operator T() const\
	{\
		init();\
		return fun_base(o1.value(0));\
	}\
};

/*! A macro to define single value function specialization that apply the function component-wise
 *
 * \param fun function name
 * \param ID function ID
 *
 */
#define CREATE_ARG_FUNC(fun_base,fun_name,OP_ID) \
\
template<typename orig,typename exp1, typename exp_2, unsigned int op1>\
__device__ __host__ inline point_expression_op<orig,point_expression_op<orig,exp1,exp_2,op1>,void, OP_ID >\
fun_name(const point_expression_op<orig,exp1,exp_2,op1> & va)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp_2,op1>,void,OP_ID> exp_sum(va);\
\
	return exp_sum;\
}\
\
\
template<unsigned int dim, typename T>\
__device__ __host__ inline point_expression_op<Point<dim,T>,Point<dim,T>,void, OP_ID >\
fun_name(const Point<dim,T> & va)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,void,OP_ID> exp_sum(va);\
\
	return exp_sum;\
}


/*! \brief Point norm operation
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 * \tparam op operation
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_NORM>
{
	const exp1 o1;

	//! Scalar value
	mutable typename std::remove_const<typename orig::coord_type>::type scal;

public:

	//! Origin type
	typedef orig orig_type;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! indicate that init must be called before value
	typedef int has_init;

	//! return type of the expression
	typedef typename orig::coord_type return_type;

	//! this operation produce a scalar as result
	static const unsigned int nvals = 1;

	//! Constructor from expression
	__device__ __host__ inline point_expression_op(const exp1 & o1)
	:o1(o1),scal(0.0)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	__device__ __host__ inline void init() const
	{
		scal = 0.0;

		for (size_t i = 0 ; i < orig::dims ; i++)
			scal += o1.value(i) * o1.value(i);

		scal = sqrt(scal);

		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > __device__ __host__ inline r_type value(size_t k) const
	{
		return scal;
	}

	template <typename T>__device__ __host__  operator T() const
	{
		init();
		return scal;
	}
};

/*! \brief Point square norm operation
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 * \tparam op operation
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_NORM2>
{
	//! expression
	const exp1 o1;

	//! result of the expression
	mutable typename std::remove_const<typename orig::coord_type>::type scal;

public:

	typedef orig orig_type;

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef typename orig::coord_type return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = 1;

	//! constructor from an expression
	__device__ __host__  inline point_expression_op(const exp1 & o1)
	:o1(o1),scal(0.0)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	__device__ __host__ inline void init() const
	{
		scal = 0.0;

		for (size_t i = 0 ; i < orig::dims ; i++)
			scal += o1.value(i) * o1.value(i);

		o1.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > __device__ __host__ inline r_type value(size_t k) const
	{
		return scal;
	}

	//! Cast operator
	template <typename T> __device__ __host__  operator T() const
	{
		init();
		return scal;
	}
};


CREATE_ARG_FUNC_CLASS(std::abs,abs,POINT_ABS)
CREATE_ARG_FUNC_CLASS(std::exp,exp,POINT_EXP)
CREATE_ARG_FUNC_CLASS(std::exp2,exp2,POINT_EXP2)
CREATE_ARG_FUNC_CLASS(std::expm1,expm1,POINT_EXPM1)
CREATE_ARG_FUNC_CLASS(std::log,log,POINT_LOG)
CREATE_ARG_FUNC_CLASS(std::log10,log10,POINT_LOG10)
CREATE_ARG_FUNC_CLASS(std::log2,log2,POINT_LOG2)
CREATE_ARG_FUNC_CLASS(std::log1p,log1p,POINT_LOG1P)
CREATE_ARG_FUNC_CLASS(std::sqrt,sqrt,POINT_SQRT)
CREATE_ARG_FUNC_CLASS(std::cbrt,cbrt,POINT_CBRT)
CREATE_ARG_FUNC_CLASS(std::sin,sin,POINT_SIN)
CREATE_ARG_FUNC_CLASS(std::cos,cos,POINT_COS)
CREATE_ARG_FUNC_CLASS(std::tan,tan,POINT_TAN)
CREATE_ARG_FUNC_CLASS(std::asin,asin,POINT_ASIN)
CREATE_ARG_FUNC_CLASS(std::acos,acos,POINT_ACOS)
CREATE_ARG_FUNC_CLASS(std::atan,atan,POINT_ATAN)
CREATE_ARG_FUNC_CLASS(std::sinh,sinh,POINT_SINH)
CREATE_ARG_FUNC_CLASS(std::cosh,cosh,POINT_COSH)
CREATE_ARG_FUNC_CLASS(std::tanh,tanh,POINT_TANH)
CREATE_ARG_FUNC_CLASS(std::asinh,asinh,POINT_ASINH)
CREATE_ARG_FUNC_CLASS(std::acosh,acosh,POINT_ACOSH)
CREATE_ARG_FUNC_CLASS(std::atanh,atanh,POINT_ATANH)
CREATE_ARG_FUNC_CLASS(std::erf,erf,POINT_ERF)
CREATE_ARG_FUNC_CLASS(std::erfc,erfc,POINT_ERFC)
CREATE_ARG_FUNC_CLASS(std::tgamma,tgamma,POINT_TGAMMA)
CREATE_ARG_FUNC_CLASS(std::lgamma,lgamma,POINT_LGAMMA)
CREATE_ARG_FUNC_CLASS(std::ceil,ceil,POINT_CEIL)
CREATE_ARG_FUNC_CLASS(std::floor,floor,POINT_FLOOR)
CREATE_ARG_FUNC_CLASS(std::trunc,trunc,POINT_TRUNC)
CREATE_ARG_FUNC_CLASS(std::round,round,POINT_ROUND)
CREATE_ARG_FUNC_CLASS(std::nearbyint,nearbyint,POINT_NEARBYINT)
CREATE_ARG_FUNC_CLASS(std::rint,rint,POINT_RINT)

/////// norm operation

/* \brief Calculate the norm of the point
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,typename exp1, typename exp2, unsigned int op1>
__device__ __host__  inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void,POINT_NORM>
norm(const point_expression_op<orig,exp1,exp2,op1> & va)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void,POINT_NORM> exp_sum(va);

	return exp_sum;
}


/* \brief norm of a scalar
 *
 * \param d scalar double or float
 *
 * \return d
 *
 */
template <typename T, typename sfinae = typename std::enable_if<std::is_fundamental<T>::value>::type >
__device__ __host__  T norm(T d)
{
	return d;
}

/* \brief Calculate the norm of the point
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename T, typename vmpl>
__device__ __host__  inline point_expression_op<Point<subar_dim<vmpl>::type::value,T>,
												point_expression<openfpm::detail::multi_array::sub_array_openfpm<T,1,vmpl>>,
												void,
												POINT_NORM>
norm(const point_expression<openfpm::detail::multi_array::sub_array_openfpm<T,1,vmpl>> & d)
{
        point_expression_op<Point<subar_dim<vmpl>::type::value,T>,
        					point_expression<openfpm::detail::multi_array::sub_array_openfpm<T,1,vmpl>>,
        					void,
        					POINT_NORM> exp_sum( (point_expression<openfpm::detail::multi_array::sub_array_openfpm<T,1,vmpl>>(d)) );

        return exp_sum;
}

/* \brief Calculate the norm of the point
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
__device__ __host__  inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,void,POINT_NORM>
norm(const point_expression<T[dim]> & d)
{
        point_expression_op<Point<dim,T>,point_expression<T[dim]>,void,POINT_NORM> exp_sum( (point_expression<T[dim]>(d)) );

        return exp_sum;
}


/* \brief Divide two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
__device__ __host__ inline point_expression_op<Point<dim,T>,Point<dim,T>,void,POINT_NORM>
norm(const Point<dim,T> & va)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,void,POINT_NORM> exp_sum(va);

	return exp_sum;
}


//////////////////////////////////////// norm2 operation //////////////////

/* \brief Calculate the norm of the point
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,typename exp1, typename exp2, unsigned int op1>
__device__ __host__ inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void,POINT_NORM2>
norm2(const point_expression_op<orig,exp1,exp2,op1> & va)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void,POINT_NORM2> exp_sum(va);

	return exp_sum;
}


/* \brief norm of a scalar
 *
 * \param d scalar double or float
 *
 * \return d
 *
 */
template <typename T>__device__ __host__ T norm2(T d)
{
	return d*d;
}


/* \brief Calculate the norm of the point
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
__device__ __host__ inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,void,POINT_NORM2>
norm2(const point_expression<T[dim]> & d)
{
        point_expression_op<Point<dim,T>,point_expression<T[dim]>,void,POINT_NORM2> exp_sum( (point_expression<T[dim]>(d)) );

        return exp_sum;
}


/* \brief Divide two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
__device__ __host__ inline point_expression_op<Point<dim,T>,Point<dim,T>,void,POINT_NORM2>
norm2(const Point<dim,T> & va)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,void,POINT_NORM2> exp_sum(va);

	return exp_sum;
}


/*! \brief General distance formula
 *
 *
 */
template <typename T, typename P> __device__ __host__  auto distance(T exp1, P exp2) -> decltype(norm(exp1 - exp2))
{
	return norm(exp1 - exp2);
}


///// Here we define point operator functions

CREATE_ARG_FUNC(std::abs,abs,POINT_ABS)
CREATE_ARG_FUNC(std::exp,exp,POINT_EXP)
CREATE_ARG_FUNC(std::exp2,exp2,POINT_EXP2)
CREATE_ARG_FUNC(std::expm1,expm1,POINT_EXPM1)
CREATE_ARG_FUNC(std::log,log,POINT_LOG)
CREATE_ARG_FUNC(std::log10,log10,POINT_LOG10)
CREATE_ARG_FUNC(std::log2,log2,POINT_LOG2)
CREATE_ARG_FUNC(std::log1p,log1p,POINT_LOG1P)
CREATE_ARG_FUNC(std::sqrt,sqrt,POINT_SQRT)
CREATE_ARG_FUNC(std::cbrt,cbrt,POINT_CBRT)
CREATE_ARG_FUNC(std::sin,sin,POINT_SIN)
CREATE_ARG_FUNC(std::cos,cos,POINT_COS)
CREATE_ARG_FUNC(std::tan,tan,POINT_TAN)
CREATE_ARG_FUNC(std::asin,asin,POINT_ASIN)
CREATE_ARG_FUNC(std::acos,acos,POINT_ACOS)
CREATE_ARG_FUNC(std::atan,atan,POINT_ATAN)
CREATE_ARG_FUNC(std::sinh,sinh,POINT_SINH)
CREATE_ARG_FUNC(std::cosh,cosh,POINT_COSH)
CREATE_ARG_FUNC(std::tanh,tanh,POINT_TANH)
CREATE_ARG_FUNC(std::asinh,asinh,POINT_ASINH)
CREATE_ARG_FUNC(std::acosh,acosh,POINT_ACOSH)
CREATE_ARG_FUNC(std::atanh,atanh,POINT_ATANH)
CREATE_ARG_FUNC(std::erf,erf,POINT_ERF)
CREATE_ARG_FUNC(std::erfc,erfc,POINT_ERFC)
CREATE_ARG_FUNC(std::tgamma,tgamma,POINT_TGAMMA)
CREATE_ARG_FUNC(std::lgamma,lgamma,POINT_LGAMMA)
CREATE_ARG_FUNC(std::ceil,ceil,POINT_CEIL)
CREATE_ARG_FUNC(std::floor,floor,POINT_FLOOR)
CREATE_ARG_FUNC(std::trunc,trunc,POINT_TRUNC)
CREATE_ARG_FUNC(std::round,round,POINT_ROUND)
CREATE_ARG_FUNC(std::nearbyint,nearbyint,POINT_NEARBYINT)
CREATE_ARG_FUNC(std::rint,rint,POINT_RINT)


#endif /* OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_FUNCTIONS_HPP_ */
