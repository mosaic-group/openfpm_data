/*
 * Point_operators.hpp
 *
 *  Created on: Jun 14, 2016
 *      Author: i-bird
 */

#ifndef OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_HPP_
#define OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_HPP_

template<unsigned int dim ,typename T> class Point;

#define POINT_SUM 1
#define POINT_SUB 2
#define POINT_MUL 3
#define POINT_DIV 4
#define POINT_MUL_POINT 5
#define POINT_NORM 6
#define POINT_NORM2 7

// Functions

#define POINT_ABS 8
#define POINT_EXP 9
#define POINT_EXP2 10
#define POINT_EXPM1 11
#define POINT_LOG 12
#define POINT_LOG10 13
#define POINT_LOG2 14
#define POINT_LOG1P 15
#define POINT_SQRT 17
#define POINT_CBRT 18
#define POINT_SIN 19
#define POINT_COS 20
#define POINT_TAN 21
#define POINT_ASIN 22
#define POINT_ACOS 23
#define POINT_ATAN 24
#define POINT_SINH 25
#define POINT_COSH 26
#define POINT_TANH 27
#define POINT_ASINH 28
#define POINT_ACOSH 29
#define POINT_ATANH 30
#define POINT_ERF 31
#define POINT_ERFC 32
#define POINT_TGAMMA 33
#define POINT_LGAMMA 34
#define POINT_CEIL 35
#define POINT_FLOOR 36
#define POINT_TRUNC 37
#define POINT_ROUND 38
#define POINT_NEARBYINT 39
#define POINT_RINT 40
#define POINT_SUB_UNI 41

constexpr unsigned int max_expr(unsigned int dim1, unsigned int dim2)
{
	return (dim1 > dim2)?dim1:dim2;
}

/*! \brief It return the dimansionality of the operation given the dimensionality of the 2 operators
 *
 * this is the defailr that return 3
 *
 * \tparam dimansionality operator1
 * \tparam dimensionality operator2
 * \tparam op operation
 *
 */
template <unsigned int op1_dim, unsigned int op2_dim, unsigned int op>
struct r_type_dim
{
	enum
	{
		value = max_expr(op1_dim,op2_dim),
	};
};

//! scalar + scalar = scalar
template <>
struct r_type_dim<1,1,POINT_SUM>
{
	enum
	{
		value = 1,
	};
};

//! scalar - scalar = scalar
template <>
struct r_type_dim<1,1,POINT_SUB>
{
	enum
	{
		value = 1,
	};
};

//! Point * Point = scalar
template <unsigned int op1_dim, unsigned int op2_dim>
struct r_type_dim<op1_dim,op2_dim,POINT_MUL_POINT>
{
	enum
	{
		value = 1,
	};
};

//! scalar * scalar = scalar
template <>
struct r_type_dim<1,1,POINT_MUL>
{
	enum
	{
		value = 1,
	};
};

//! scalar / scalar = scalar
template <>
struct r_type_dim<1,1,POINT_DIV>
{
	enum
	{
		value = 1,
	};
};

/*! \brief Return type of the expression
 *
 * \tparam r dimension of the return type
 * \tparam orig original type
 *
 */
template <unsigned int r, typename orig>
struct r_type_p
{
	typedef orig type;
};

/*! \brief Return type of the expression
 *
 * \tparam orig original type
 *
 */
template <typename orig>
struct r_type_p<1,orig>
{
	typedef typename orig::coord_type type;
};


/*! \brief Main class that encapsulate a constant number used in a point expression
 *
 *
 */
template<typename T>
class point_expression
{
	//! constant
	T d;

public:

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = 1;

	//! constructor from a value
	inline point_expression(T & d)
	:d(d)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k coordinate to evaluate
	 *
	 * \return It just return the velue set in the constructor
	 *
	 */
	inline T value(const size_t k) const
	{
		return d;
	}
};



/*! \brief Unknown operation specialization
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 * \tparam op operation
 *
 */
template <typename orig, typename exp1, typename exp2, unsigned int op>
class point_expression_op
{

};

/*! \brief Sum operation
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_SUM>
{
	//! first expression
	const exp1 o1;
	//! second expression
	const exp2 o2;

public:

	//! original type of the point expression
	typedef orig orig_type;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! indicate that init must be called before value
	typedef int has_init;

	//! return type of the expression
	typedef typename r_type_p<r_type_dim<exp1::nvals,exp2::nvals,POINT_SUM>::value,orig >::type return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_SUM>::value;

	//! Constructor from 2 point expressions
	inline point_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression at coordinate k
	 *
	 * \param k coordinate
	 *
	 * \return the value of the expression for the coordinate k
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return o1.value(k) + o2.value(k);
	}

	/*! \brief conversion of the class to double or float or ...
	 *
	 *
	 */
	template<typename T, typename test=typename boost::disable_if_c< std::is_same<T,orig>::value || exp1::nvals != 1 || exp2::nvals != 1 >::type  > inline operator T() const
	{
		init();
		return o1.value(0) + o2.value(0);
	}
};

/*! \brief Subtraction operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename orig,typename exp1, typename exp2>
class point_expression_op<orig, exp1,exp2,POINT_SUB>
{
	//! expression 1
	const exp1 o1;
	//! expression 2
	const exp2 o2;

public:

	//! Original type
	typedef orig orig_type;

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef orig return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_SUB>::value;

	//! constructor from 2 expressions
	inline point_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression at coordinate k
	 *
	 * \param k coordinate
	 *
	 * \return the evaluate expression at coordinate k
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return o1.value(k) - o2.value(k);
	}

	/*! \brief conversion of the class to double or float or ...
	 *
	 *
	 */
	template<typename T, typename test=typename boost::disable_if_c< std::is_same<T,orig>::value || exp1::nvals != 1 || exp2::nvals != 1 >::type  > operator T() const
	{
		init();
		return o1.value(0) - o2.value(0);
	}
};

template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2, POINT_SUB_UNI >
{
	//! expression
	const exp1 o1;

	//! scalar value produced by the expression
	mutable typename orig::coord_type scal;

public:

	//! original type
	typedef orig orig_type;

	//! indicate that is an expression
	typedef int is_expression;

	//! indicate that this class has an init function
	typedef int has_init;

	//! return type of the expression evaluation
	typedef typename orig::coord_type return_type;

	//! result dimensionality of this expression
	static const unsigned int nvals = exp1::nvals;

	//! constructor from expression
	inline point_expression_op(const exp1 & o1)
	:o1(o1),scal(0.0)
	{}

	//! initialize the the expression
	inline void init() const
	{
		o1.init();
	}

	//! evaluate the expression
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return -(o1.value(k));
	}

	//! casting to a type T
	template <typename T, typename check = typename std::enable_if<!std::is_same<T,orig>::value >::type  >operator T() const
	{
		init();
		return -(o1.value(0));
	}
};


/*! \brief Multiplication operation
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 * \tparam op operation
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_MUL_POINT>
{
	//! first expression
	const exp1 o1;
	//! second expression
	const exp2 o2;

	mutable typename orig::coord_type scal;

public:

	typedef orig orig_type;

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef typename orig::coord_type return_type;

	//! this operation produce a scalar as result
	static const unsigned int nvals = 1;

	//! constructor from 2 expressions
	inline point_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2),scal(0.0)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
		for (size_t i = 0 ; i < orig::dims ; i++)
			scal += o1.value(i) * o2.value(i);

		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return scal;
	}

	//! cast to other type
	template<typename T, typename test=typename boost::disable_if_c< std::is_same<T,orig>::value >::type > operator T() const
	{
		init();
		return scal;
	}
};


/*! \brief Multiplication operation
 *
 * \tparam orig original type
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_MUL>
{
	//! expression 1
	const exp1 o1;
	//! expression 2
	const exp2 o2;

public:

	//! origin type
	typedef orig orig_type;

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef orig return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_MUL>::value;

	//! constructor from 2 expression
	inline point_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param key where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return o1.value(k) * o2.value(k);
	}

	/*! \brief conversion of the class to double or float or ...
	 *
	 *
	 */
	template<typename T, typename test=typename boost::disable_if_c< std::is_same<T,orig>::value || exp1::nvals != 1 || exp2::nvals != 1 >::type  > operator T() const
	{
		init();
		return o1.value(0) * o2.value(0);
	}

};

/*! \brief Division operation
 *
 * \tparam exp1 expression1
 * \tparam exp2 expression2
 *
 */
template <typename orig, typename exp1, typename exp2>
class point_expression_op<orig,exp1,exp2,POINT_DIV>
{
	//! expression 1
	const exp1 o1;
	//! expression 2
	const exp2 o2;

public:

	//! original type
	typedef orig orig_type;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! indicate that init must be called before value
	typedef int has_init;

	//! return type of the expression
	typedef orig return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_DIV>::value;

	//! constructor from expression 1 and expression 2
	inline point_expression_op(const exp1 & o1, const exp2 & o2)
	:o1(o1),o2(o2)
	{}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
		o1.init();
		o2.init();
	}

	/*! \brief Evaluate the expression
	 *
	 * \param k where to evaluate the expression
	 *
	 */
	template<typename r_type=typename std::remove_reference<decltype(o1.value(0))>::type > inline r_type value(size_t k) const
	{
		return o1.value(k) / o2.value(k);
	}


	/*! \brief conversion of the class to double or float or ...
	 *
	 *
	 */
	template<typename T, typename test=typename boost::disable_if_c< std::is_same<T,orig>::value || exp1::nvals != 1 || exp2::nvals != 1 >::type  > operator T() const
	{
		init();
		return o1.value(0) / o2.value(0);
	}
};

/*! \brief Transform an array into a point expression
 *
 * \param array
 *
 * \return an object that can be used into an expression
 *
 */
template<unsigned int dim, typename T> point_expression<T[dim]> getExpr(T (& a)[dim])
{
	return point_expression<T[dim]>(a);
}

/*! \brief MACRO that define operator for point template expression parsing
 *
 *
 */

#define CREATE_POINT_OPERATOR(operator_name,OP_ID) \
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T[dim]>,OP_ID>\
operator_name(const Point<dim,T> & va, const point_expression<T[(unsigned int)dim]> & vb)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T[dim]>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,Point<dim,T>,OP_ID>\
operator_name(const point_expression<T[(unsigned int)dim]> & va, const Point<dim,T> & vb)\
{\
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,Point<dim,T>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<double>,OP_ID>\
operator_name(const point_expression<T[dim]> & va, double d)\
{\
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<double>,OP_ID> exp_sum(va,point_expression<double>(d));\
\
	return exp_sum;\
}\
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,point_expression<double>,point_expression<T[dim]>,OP_ID>\
operator_name(double d, const point_expression<T[dim]> & va)\
{\
	point_expression_op<Point<dim,T>,point_expression<double>,point_expression<T[dim]>,OP_ID> exp_sum(point_expression<double>(d),va);\
\
	return exp_sum;\
}\
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<T[dim]>,OP_ID>\
operator_name(const point_expression<T[dim]> & va, const point_expression<T[dim]> & vb)\
{\
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<T[dim]>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T[dim]>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression<T[dim]> & vb)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T[dim]>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>\
inline point_expression_op<orig,point_expression<T[dim]>,point_expression_op<orig,exp1,exp2,op1>,OP_ID>\
operator_name(const point_expression<T[dim]> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)\
{\
	point_expression_op<orig,point_expression<T[dim]>,point_expression_op<orig,exp1,exp2,op1>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
\
template<unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,OP_ID>\
operator_name(const Point<dim,T> & va, const Point<dim,T> & vb)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename orig,typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>\
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,OP_ID>\
operator_name(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T>\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, T d)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,OP_ID> exp_sum(va,d);\
\
	return exp_sum;\
}\
\
template<typename orig,typename exp1 , typename exp2, unsigned int op1, typename T>\
inline point_expression_op<orig,point_expression<T>,point_expression_op<orig,exp1,exp2,op1>,OP_ID>\
operator_name(T d, const point_expression_op<orig,exp1,exp2,op1> & vb)\
{\
	point_expression_op<orig,point_expression<T>,point_expression_op<orig,exp1,exp2,op1>,OP_ID> exp_sum(d,vb);\
\
	return exp_sum;\
}\
\
template<typename orig,typename exp1 , typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,OP_ID> exp_sum(va,vb);\
\
	return exp_sum;\
}\
\
template<unsigned int dim , typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >\
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,OP_ID>\
operator_name(const Point<dim,T> & va, T d)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,OP_ID> exp_sum(va,point_expression<T>(d));\
\
	return exp_sum;\
}\
\
template<unsigned int dim , typename T>\
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,OP_ID>\
operator_name(const Point<dim,T> & va, double d)\
{\
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,OP_ID> exp_sum(va,point_expression<double>(d));\
\
	return exp_sum;\
}\
\
template<unsigned int dim , typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >\
inline point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,OP_ID>\
operator_name(T d, const Point<dim,T> & vb)\
{\
	point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,OP_ID> exp_sum(point_expression<T>(d),vb);\
\
	return exp_sum;\
}\
\
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type>\
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,Point<dim,T>,OP_ID>\
operator_name(T (& d)[dim], const Point<dim,T> & vb)\
{\
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,Point<dim,T>,OP_ID> exp_sum(point_expression<T[dim]>(d),vb);\
\
	return exp_sum;\
}\
\
template<unsigned int dim , typename T>\
inline point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,OP_ID>\
operator_name(double d, const Point<dim,T> & vb)\
{\
	point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,OP_ID> exp_sum(point_expression<double>(d),vb);\
\
	return exp_sum;\
}\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, T d)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,OP_ID> exp_sum(va,point_expression<T>(d));\
\
	return exp_sum;\
}\
\
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T>\
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,OP_ID>\
operator_name(const point_expression_op<orig,exp1,exp2,op1> & va, double d)\
{\
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,OP_ID> exp_sum(va,point_expression<double>(d));\
\
	return exp_sum;\
}


CREATE_POINT_OPERATOR(operator+,POINT_SUM)
CREATE_POINT_OPERATOR(operator-,POINT_SUB)
CREATE_POINT_OPERATOR(operator/,POINT_DIV)

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< ! std::is_same<T,orig>::value >::type >
inline T &
operator+=(T & d, const point_expression_op<orig,exp1,exp2,op1> & va)
{
	va.init();
	d += va.value(0);

	return d;
}



/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T>
inline T &
operator+=(orig & d, const point_expression_op<orig,exp1,exp2,op1> & va)
{
	va.init();

	d = d + va;

	return d;
}

/* \brief minus points expression
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void, POINT_SUB_UNI >
operator-(const point_expression_op<orig,exp1,exp2,op1> & va)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,void,POINT_SUB_UNI> exp_sum(va);

	return exp_sum;
}

/* \brief minus points expression
 *
 * \param va point expression one
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,void, POINT_SUB_UNI >
operator-(const Point<dim,T> & va)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,void,POINT_SUB_UNI> exp_sum(va);

	return exp_sum;
}


///////////////////// DEFINITION OF THE OPERATOR* ////////////////
///////////////// Template expression parsing ////////////////////


template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<double>,POINT_MUL>
operator*(const point_expression<T[dim]> & va, double d)
{
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<double>,POINT_MUL> exp_sum(va,point_expression<double>(d));

	return exp_sum;
}

template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,point_expression<double>,point_expression<T[dim]>,POINT_MUL>
operator*(double d, const point_expression<T[dim]> & va)
{
	point_expression_op<Point<dim,T>,point_expression<double>,point_expression<T[dim]>,POINT_MUL> exp_sum(point_expression<double>(d),va);

	return exp_sum;
}

template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<T[dim]>,POINT_MUL_POINT>
operator*(const point_expression<T[dim]> & va, const point_expression<T[dim]> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<T[dim]>,point_expression<T[dim]>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type>
inline point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_MUL>
operator*(T d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_MUL> exp_sum(point_expression<T>(d),vb);

	return exp_sum;
}

template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T[dim]>,POINT_MUL_POINT>
operator*(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression<T[dim]> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T[dim]>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<orig,point_expression<T[dim]>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL_POINT>
operator*(const point_expression<T[dim]> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,point_expression<T[dim]>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_MUL>
operator*(double d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_MUL> exp_sum(point_expression<double>(d),vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_MUL>
operator*(const Point<dim,T> & va, T d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_MUL> exp_sum(va,point_expression<T>(d));

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_MUL>
operator*(const Point<dim,T> & va, double d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_MUL> exp_sum(va,point_expression<double>(d));

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_MUL_POINT>
operator*(const Point<dim,T> & va, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, unsigned int dim, typename T, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL_POINT>
operator*(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,unsigned int dim, typename T, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_MUL_POINT>
operator*(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply a point expression with a number
 *
 * \param d number
 * \param vb point expression
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename T, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,point_expression<T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL>
operator*(T d, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,point_expression<T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL> exp_sum(point_expression<T>(d),vb);

	return exp_sum;
}


/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param d constant value
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename T>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_MUL>
operator*(const point_expression_op<orig,exp1,exp2,op1> & va, T d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_MUL> exp_sum(va,point_expression<T>(d));

	return exp_sum;
}


/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_MUL_POINT>
operator*(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_MUL_POINT> exp_sum(va,vb);

	return exp_sum;
}

////////////////////////////// Point wise multiplication ///////////////////////

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_MUL>
pmul(const Point<dim,T> & va, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, unsigned int dim, typename T, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL>
pmul(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,unsigned int dim, typename T, typename exp1, typename exp2, unsigned int op1>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_MUL>
pmul(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_MUL> exp_sum(va,vb);

	return exp_sum;
}

/* \brief Multiply two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_MUL>
pmul(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_MUL> exp_sum(va,vb);

	return exp_sum;
}


/*! \brief Specialization for an array of dimension dim
 *
 * \tparam type
 *
 */
template<typename T, unsigned int dim>
class point_expression<T[dim]>
{
	T (& d)[dim];

public:

	// indicate that init must be called before value
	typedef int has_init;

	// indicate that this class encapsulate an expression
	typedef int is_expression;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = dim;

	inline point_expression(T (& d)[dim])
	:d(d)
	{
	}

	/*! \brief Operator= for point expression
	 *
	 * \tparam orig origin type
	 * \tparam exp1 expression 1
	 * \tparam exp2 expression 2
	 * \tparam op operation
	 *
	 * \param point expression
	 *
	 * \return a point expression
	 *
	 */
	template<typename orig, typename exp1, typename exp2, unsigned int op> point_expression<T[dim]> & operator=(const point_expression_op<orig,exp1,exp2,op> & p_exp)
	{
		p_exp.init();

		for (size_t i = 0; i < dim ; i++)
			d[i] = p_exp.value(i);

		return *this;
	}

	/*! \brief This function must be called before value
	 *
	 * it calculate the scalar product before return the values
	 *
	 */
	inline void init() const
	{
	}

	/*! \brief Evaluate the expression at coordinate k
	 *
	 * It just return the value set in the constructor
	 *
	 * \param k coordinate
	 *
	 * \return the value
	 *
	 */
	inline T value(const size_t k) const
	{
		return d[k];
	}
};


#include "Point_operators_functions.hpp"

#endif /* OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_HPP_ */
