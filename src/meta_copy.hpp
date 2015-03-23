#ifndef META_COPY_HPP
#define META_COPY_HPP

/*! \brief This class is an helper to copy scalar and compile-time array element
 *
 * This class is an helper to copy scalar and compile-time array element
 *
 * Usage:
 *
 * float src[3];
 * float dst[3];
 *
 * meta_copy<T[3]> cp(src,dst);
 *
 *
 */

/*! \brief Scalar copy
 *
 * Scalar copy
 *
 * \param T type of element to copy
 *
 * This specialization pattern is chosen in case of scalar
 *
 */
template<typename T>
struct meta_copy
{
	meta_copy(const T & src, T & dst)
	{
		dst = src;
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct meta_copy<T[N1]>
{
	meta_copy(const T src[N1], T dst[N1])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			dst[i1] = src[i1];
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct meta_copy<T[N1][N2]>
{
	meta_copy(const T src[N1][N2], T dst[N1][N2])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				dst[i1][i2] = src[i1][i2];
			}
		}
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct meta_copy<T[N1][N2][N3]>
{
	meta_copy(const T src[N1][N2][N3], T dst[N1][N2][N3])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					dst[i1][i2][i3] = src[i1][i2][i3];
				}
			}
		}
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct meta_copy<T[N1][N2][N3][N4]>
{
	meta_copy(const T src[N1][N2][N3][N4], T dst[N1][N2][N3][N4])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						dst[i1][i2][i3][i4] = src[i1][i2][i3][i4];
					}
				}
			}
		}
	}
};

//! Partial specialization for N=5
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5>
struct meta_copy<T[N1][N2][N3][N4][N5]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5], T dst[N1][N2][N3][N4][N5])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							dst[i1][i2][i3][i4][i5] = src[i1][i2][i3][i4][i5];
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=6
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6>
struct meta_copy<T[N1][N2][N3][N4][N5][N6]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5][N6], T dst[N1][N2][N3][N4][N5][N6])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								dst[i1][i2][i3][i4][i5][i6] = src[i1][i2][i3][i4][i5][i6];
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=7
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5][N6][N7], T dst[N1][N2][N3][N4][N5][N6][N7])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									dst[i1][i2][i3][i4][i5][i6][i7] = src[i1][i2][i3][i4][i5][i6][i7];
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=8
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5][N6][N7][N8], T dst[N1][N2][N3][N4][N5][N6][N7][N8])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										dst[i1][i2][i3][i4][i5][i6][i7][i8] = src[i1][i2][i3][i4][i5][i6][i7][i8];
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=9
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8][N9]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5][N6][N7][N8][N9], T dst[N1][N2][N3][N4][N5][N6][N7][N8][N9])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										for (size_t i9 = 0 ; i9 < N9 ; i9++)
										{
											dst[i1][i2][i3][i4][i5][i6][i7][i8][i9] = src[i1][i2][i3][i4][i5][i6][i7][i8][i9];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

//! Partial specialization for N=10
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9, size_t N10>
struct meta_copy<T[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10]>
{
	meta_copy(const T src[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10], T dst[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10])
	{
		for (size_t i1 = 0 ; i1 < N1 ; i1++)
		{
			for (size_t i2 = 0 ; i2 < N2 ; i2++)
			{
				for (size_t i3 = 0 ; i3 < N3 ; i3++)
				{
					for (size_t i4 = 0 ; i4 < N4 ; i4++)
					{
						for (size_t i5 = 0 ; i5 < N5 ; i5++)
						{
							for (size_t i6 = 0 ; i6 < N6 ; i6++)
							{
								for (size_t i7 = 0 ; i7 < N7 ; i7++)
								{
									for (size_t i8 = 0 ; i8 < N8 ; i8++)
									{
										for (size_t i9 = 0 ; i9 < N9 ; i9++)
										{
											for (size_t i10 = 0 ; i10 < N10 ; i10++)
											{
												dst[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10] = src[i1][i2][i3][i4][i5][i6][i7][i8][i9][i10];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
};

#endif
