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

/*! \brief It return the dimansionality of the operation given the dimensionality of the 2 operators
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
		value = 3,
	};
};

template <>
struct r_type_dim<1,1,POINT_SUM>
{
	enum
	{
		value = 1,
	};
};

template <>
struct r_type_dim<1,1,POINT_SUB>
{
	enum
	{
		value = 1,
	};
};

template <unsigned int op1_dim, unsigned int op2_dim>
struct r_type_dim<op1_dim,op2_dim,POINT_MUL_POINT>
{
	enum
	{
		value = 1,
	};
};

template <>
struct r_type_dim<1,1,POINT_MUL>
{
	enum
	{
		value = 1,
	};
};

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


/*! \brief Main class that encapsulate a constant number
 *
 * \param prp no meaning
 *
 */
template<typename T>
class point_expression
{
	T d;

public:

	// indicate that init must be called before value
	typedef int has_init;

	// indicate that this class encapsulate an expression
	typedef int is_expression;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = 1;

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
	 * It just return the velue set in the constructor
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
	const exp1 o1;
	const exp2 o2;

public:

	typedef orig orig_type;

	// indicate that this class encapsulate an expression
	typedef int is_expression;

	// indicate that init must be called before value
	typedef int has_init;

	//! return type of the expression
	typedef typename r_type_p<r_type_dim<exp1::nvals,exp2::nvals,POINT_SUM>::value,orig >::type return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_SUM>::value;

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
	 * \param k dimension
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
	const exp1 o1;
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
	 * \param k dimension
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
	const exp1 o1;
	const exp2 o2;

	mutable typename orig::coord_type scal;

public:

	typedef orig orig_type;

	// indicate that init must be called before value
	typedef int has_init;

	// indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef typename orig::coord_type return_type;

	//! this operation produce a scalar as result
	static const unsigned int nvals = 1;

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
		scal = 0.0;

		for (size_t i = 0 ; i < orig::dims ; i++)
			scal += o1.value(i) * o2.value(i);

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
		return scal;
	}

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
	const exp1 o1;
	const exp2 o2;

public:

	typedef orig orig_type;

	//! indicate that init must be called before value
	typedef int has_init;

	//! indicate that this class encapsulate an expression
	typedef int is_expression;

	//! return type of the expression
	typedef orig return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_MUL>::value;

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
	const exp1 o1;
	const exp2 o2;

public:

	typedef orig orig_type;

	// indicate that this class encapsulate an expression
	typedef int is_expression;

	// indicate that init must be called before value
	typedef int has_init;

	//! return type of the expression
	typedef orig return_type;

	//! this operation produce a vector as result of size dims
	static const unsigned int nvals = r_type_dim<exp1::nvals,exp2::nvals,POINT_DIV>::value;

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


/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_SUM>
operator+(const Point<dim,T> & va, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_SUM>
operator+(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,typename exp1 , typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_SUM>
operator+(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,typename exp1 , typename exp2, unsigned int op1, typename exp3 , typename exp4, unsigned int op2>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_SUM>
operator+(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_SUM> exp_sum(va,vb);

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim , typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_SUM>
operator+(const Point<dim,T> & va, T d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_SUM> exp_sum(va,point_expression<T>(d));

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim , typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_SUM>
operator+(const Point<dim,T> & va, double d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_SUM> exp_sum(va,point_expression<double>(d));

	return exp_sum;
}

/* \brief sum two point expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim , typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_SUM>
operator+(T d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_SUM> exp_sum(point_expression<T>(d),vb);

	return exp_sum;
}



/* \brief sum two point expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim , typename T>
inline point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_SUM>
operator+(double d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_SUM> exp_sum(point_expression<double>(d),vb);

	return exp_sum;
}

/* \brief sum two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1 , typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_SUM>
operator+(const point_expression_op<orig,exp1,exp2,op1> & va, T d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_SUM> exp_sum(va,point_expression<T>(d));

	return exp_sum;
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
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_SUM>
operator+(const point_expression_op<orig,exp1,exp2,op1> & va, double d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_SUM> exp_sum(va,point_expression<double>(d));

	return exp_sum;
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
operator+=(T & d, const point_expression_op<orig,exp1,exp2,op1> & va)
{
	va.init();
	d += va.value(0);

	return d;
}

/* \brief subtract two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_SUB>
operator-(const Point<dim,T> & va, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_SUB> exp_sum(va,vb);

	return exp_sum;
}


/* \brief subtract two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_SUB>
operator-(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, unsigned int dim, typename T>
inline point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_SUB>
operator-(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,Point<dim,T>, point_expression_op<orig,exp1,exp2,op1>,POINT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_SUB>
operator-(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_SUB> exp_sum(va,vb);

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_SUB>
operator-(const Point<dim,T> & va, T d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_SUB> exp_sum(va,point_expression<T>(d));

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_SUB>
operator-(const Point<dim,T> & va, double d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_SUB> exp_sum(va,point_expression<double>(d));

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type>
inline point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_SUB>
operator-(T d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_SUB> exp_sum(point_expression<T>(d),vb);

	return exp_sum;
}

/* \brief subtract two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<unsigned int dim, typename T>
inline point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_SUB>
operator-(double d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_SUB> exp_sum(point_expression<double>(d),vb);

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
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type>
inline point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_MUL>
operator*(T d, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,point_expression<T>,Point<dim,T>,POINT_MUL> exp_sum(point_expression<T>(d),vb);

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
 * \param d constant value
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_MUL>
operator*(const point_expression_op<orig,exp1,exp2,op1> & va, double d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_MUL> exp_sum(va,point_expression<double>(d));

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

/////////////////////////////////////////////////

/* \brief Divide two points expression
 *
 * \param va point expression one
 * \param vb point expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig,typename exp1, typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_DIV>
operator/(const point_expression_op<orig,exp1,exp2,op1> & va, T d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_DIV> exp_sum(va,point_expression<T>(d));

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
template<typename orig,typename exp1, typename exp2, unsigned int op1, typename T>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_DIV>
operator/(const point_expression_op<orig,exp1,exp2,op1> & va, double d)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<double>,POINT_DIV> exp_sum(va,point_expression<double>(d));

	return exp_sum;
}

/* \brief Divide two distributed vector expression
 *
 * \param va vector expression one
 * \param vb vector expression two
 *
 * \return an object that encapsulate the expression
 *
 */
template<typename orig, typename exp1, typename exp2, unsigned int op1, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type >
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression<T>,POINT_DIV>
operator/(T d, const point_expression_op<orig,exp1,exp2,op1> & va)
{
	point_expression_op<orig,point_expression<T>,point_expression_op<orig,exp1,exp2,op1>,POINT_DIV> exp_sum(point_expression<T>(d),va);

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
inline point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_DIV>
operator/(double d, const Point<dim,T> & va)
{
	point_expression_op<Point<dim,T>,point_expression<double>,Point<dim,T>,POINT_DIV> exp_sum(point_expression<double>(d),va);

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
template<unsigned int dim, typename T, typename check=typename std::enable_if< !std::is_same<T,double>::value >::type>
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_DIV>
operator/(const Point<dim,T> & va, T d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<T>,POINT_DIV> exp_sum(va,point_expression<T>(d));

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
inline point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_DIV>
operator/(const Point<dim,T> & va, double d)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,point_expression<double>,POINT_DIV> exp_sum(va,point_expression<double>(d));

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
inline point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_DIV>
operator/(const Point<dim,T> & va, const Point<dim,T> & vb)
{
	point_expression_op<Point<dim,T>,Point<dim,T>,Point<dim,T>,POINT_DIV> exp_sum(va,vb);

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
template<typename orig,unsigned int dim, typename T, typename exp1,typename exp2, unsigned int op1>
inline point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_DIV>
operator/(const Point<dim,T> & va, const point_expression_op<orig,exp1,exp2,op1> & vb)
{
	point_expression_op<orig,Point<dim,T>,point_expression_op<orig,exp1,exp2,op1>,POINT_DIV> exp_sum(va,vb);

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
template<typename orig,unsigned int dim, typename T, typename exp1,typename exp2, unsigned int op1>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_DIV>
operator/(const point_expression_op<orig,exp1,exp2,op1> & va, const Point<dim,T> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,Point<dim,T>,POINT_DIV> exp_sum(va,vb);

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
template<typename orig,typename exp1,typename exp2, unsigned int op1, typename exp3, typename exp4, unsigned int op2>
inline point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_DIV>
operator/(const point_expression_op<orig,exp1,exp2,op1> & va, const point_expression_op<orig,exp3,exp4,op2> & vb)
{
	point_expression_op<orig,point_expression_op<orig,exp1,exp2,op1>,point_expression_op<orig,exp3,exp4,op2>,POINT_DIV> exp_sum(va,vb);

	return exp_sum;
}


#include "Point_operators_functions.hpp"

#endif /* OPENFPM_DATA_SRC_SPACE_SHAPE_POINT_OPERATORS_HPP_ */
