#ifndef HYPERCUBE_HPP
#define HYPERCUBE_HPP

#include "mathutil.hpp"
#include "Grid/grid_sm.hpp"
#include "Grid/comb.hpp"
#include "mathutil.hpp"

template<unsigned int dim, unsigned int subdim> class SubHyperCube;

/*! \brief output stream overload for printing
 *
 * \param str ostream
 * \param c combination to print
 *
 */
template<unsigned int dim> std::ostream& operator<<(std::ostream& str, const comb<dim> & c)
{
	// print the combination of ostream
	for (size_t i = 0 ; i < dim-1 ; i++)
		str <<  c.c[i] << ";";

	str << c.c[dim-1];

    return str;
}

/*! \brief This class calculate elements of the hyper-cube
 *
 * This class give you a set of utility functions for the hyper-cube like getting
 * number of faces, number of edge, number of vertex, or in general number of elements
 * of dimension d, position of each element
 *
 * * 0d Hyper-cube vertex
 * * 1d Hypercube segment
 * * 2d Hypercube square
 * * 3d Hypercube Cube
 * ...
 *
 * \tparam dim dimensionality of the Hyper-cube
 *
 * ### Get vertex and edge on a line
 * \snippet hypercube_unit_test.hpp Get vertex and edge on a line
 * ### Get vertex edge and surfaces of a square
 * \snippet hypercube_unit_test.hpp  Get vertex edge and surfaces of a square
 * ### Get vertex edge surfaces and volumes of a cube
 * \snippet hypercube_unit_test.hpp Get vertex edge surfaces and volumes of a cube
 *
 */

template<unsigned int dim>
class HyperCube
{
public:

	/*! \brief Get the number of Elements of dimension d
	 *
	 * \param d dimensionality of the element
	 * \return the number of elements of dimension d
	 *
	 */
	static size_t getNumberOfElements_R(size_t d)
	{
		// The formula to calculate the number of element of size d is given by
		//
		// 2^(dim-d) * C(dim,d) with C(dim,d) = dim!/(d!(dim-d)!)

		return pow(2,dim-d) * openfpm::math::C(dim,d);
	}

	/*! \brief Get the sum of the number of elements from d to d_t (included)
	 *
	 * \param d_t
	 * \return the sum of the number of elements from d to d_t
	 *
	 */

	static size_t getNumberOfElementsTo_R(size_t d_t)
	{
		size_t N_ele = 0;

		for (size_t i = 0 ; i <= d_t ; i++)
			N_ele += getNumberOfElements_R(i);

		return N_ele;
	}

	/*! brief Calculate the position (combinations) of all the elements of size d
	 *
	 * \param d dimensionality of the object position
	 * \return all the combinations
	 *
	 */

	static std::vector<comb<dim>> getCombinations_R(size_t d)
	{
		// Create an Iterator_g_const
		// And a vector that store all the combination

		std::vector<comb<dim>> v;
		Iterator_g_const it(dim-d,dim);

		// for each non-zero elements
		// basically the key will store the position of the
		// non zero elements, while BinPermutations will
		// fill the array of all the permutations
		//

		while (it.isNext())
		{
			grid_key_dx_r & key = it.get();

			// Calculate the permutation

			BinPermutations(key,v);
			++it;
		}

		// case when d == dim

		if (d == dim)
		{
			comb<dim> c;
			c.zero();

			v.push_back(c);
		}

		// return the combinations

		return v;
	}

	/*! \brief Binary permutations
	 *
	 * Fill v with all the possible binary permutations
	 * it produce 2^(pos.getDim()) Permutations
	 *
	 * Example
	 *
	 * if getDim() is 2
	 *
	 * it produce 4 configuration
	 *
	 * (1,1) (1,-1) (-1,1) (-1,-1)
	 *
	 * and fill the number in the position indicated by Iterator_g_const
	 *
	 * \param pos slots inside comb to fill with all permutations
	 * \param v vector to fill with the permutations
	 *
	 */

	static void BinPermutations(grid_key_dx_r & pos, std::vector<comb<dim>> & v)
	{
		size_t p_stop = pow(2,pos.getDim());
		for (size_t p = 0 ; p < p_stop ; p++)
		{
			// Create a new permutation
			struct comb<dim> c;

			// zero the combination
			c.zero();

			// the bit of p give the permutations 0 mean bottom 1 mean up
			// Fill the combination based on the bit of p

			for (size_t k = 0 ; k < pos.getDim() ; k++)
			{
				// bit of p
				bool bit = (p >> k) & 0x01;

				// Fill the combination
				c.c[pos.get(k)] = (bit == 0)? 1 : -1;
			}

			// save the combination

			v.push_back(c);
		}
	}

	static SubHyperCube<dim,dim-1> getSubHyperCube(int d)
	{
		SubHyperCube<dim,dim-1> sub;

		return sub;
	}

//	std::vector<comb<dim>> gm;

	/** \brief Linearize the combination
	 *
	 *	\param c given combination
	 *  \return the linearized combination
	 *
	 */

	static size_t LinId(comb<dim> & c)
	{
		// calculate the numbers of non-zero
		size_t d = 0;

		size_t pos_n_zero[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c.c[i] != 0)
			{d++;}
		}

		// Get the position of the non-zero
		size_t pn_zero = 0;
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (c.c[i] != 0)
			{
				pos_n_zero[d-pn_zero-1] = i;
				pn_zero++;
			}
		}

		// The combination linearization is given by
		size_t lin_id = 0;

		// Cumulative value
		long int val = 0;
		long int cum_val = 0;
		for(long int i = d - 1; i >= 0 ; i--)
		{
			// check the out-of-bound, outside is assumed to be -1 so  (- pos_n_zero[i+1] - 1) = 0
			if (i+1 < (long int)d)
				val = pos_n_zero[i] - pos_n_zero[i+1] - 1;
			else
				val = pos_n_zero[i];

			for (long int j = 0 ; j < (long int)val; j++)
			{
				// C is not safe check the limit
				if (((long int)dim)-cum_val-j-1 >= 0 && i > 0 && ((long int)dim)-cum_val-j >= i)
					lin_id += openfpm::math::C(dim-cum_val-j-1,i);
				else
					lin_id += 1;
			}

			cum_val += (val + 1);
		}

		// multiply for the permutation
		lin_id *= pow(2,d);

		// calculate the permutation position
		size_t id = 0;

		for (size_t i = 0 ; i < d ; i++)
		{
			if (c.c[pos_n_zero[i]] == -1)
			{id = id | (1 << i);}
		}

		// return the correct id

		return lin_id + id;
	}

	/** \brief isPositive return if the combination d is a positive or a negative
	 *
	 * For an hyper-cube of dimension dim we have 2*dim combinations half on positive direction
	 * half on negative direction, the function check
	 * if the d combination is negative or positive
	 *
	 * \param d
	 *
	 */
	static bool isPositive(size_t d)
	{
		return (d % 2) == 0;
	}

	/*! \brief return the combination of the positive face on direction d
	 *
	 * \param d direction
	 *
	 * \return id of the combination
	 *
	 */
	static int positiveFace(int d)
	{
		return d * 2;
	}

	/*! \brief Return the combination of the negative face on direction d
	 *
	 * \param d direction
	 *
	 * \return id of the combination
	 *
	 */
	static int negativeFace(int d)
	{
		return d * 2 + 1;
	}
};

/*! \brief This represent a sub-hyper-cube of an hyper-cube like a face or an edge of a cube
 *
 * It give a set of utility function to work with sub-hyper-cubes like the hyper-cube
 *
 * \tparam dimensionality of the hyper-cube
 * \tparam dimensionality of the sub-hyper-cube
 *
 * ### Getting the surfaces of the cube
 * \snippet hypercube_unit_test.hpp Getting the surfaces of the cube
 * ### Getting the vertices of the surfaces of the cube
 * \snippet hypercube_unit_test.hpp Getting the vertices of the surfaces of the cube
 * ### Getting the edges of the surfaces of the cube
 * \snippet hypercube_unit_test.hpp Getting the edges of the surfaces of the cube
 *
 */

template<unsigned int dim, unsigned int subdim>
class SubHyperCube : public HyperCube<subdim>
{
public:

	/*! brief Calculate the position (combinations) of all the elements of size d in the sub-hyper-cube
	 *
	 * \param c identify the position of the sub-hyper-cube in the hypercube
	 * \param d dimensionality of the objects
	 * \return all the combinations
	 *
	 */
	static std::vector<comb<dim>> getCombinations_R(comb<dim> c, int d)
	{
#ifdef DEBUG
		if (c.n_zero() < d)
		{
			std::cerr << "Error SubHyperCube: " << __FILE__ << " " << __LINE__ << " the hyper-cube selected must have dimensionality bigger than the dimensionality of the requested combinations, or the number of zero in c must me bigger than d" << "\n";
		}
#endif

		// if sub-dim == 0 return c

		if (subdim == 0)
		{
			std::vector<comb<dim>> vc;
			vc.push_back(c);

			return vc;
		}

		// Create an Iterator_g_const
		// And a vector that store all the combination

		std::vector<comb<subdim>> v = HyperCube<subdim>::getCombinations_R(d);

		// Create new combinations
		std::vector<comb<dim>> vc(v.size());

		// for each combination
		for (size_t i = 0 ; i < v.size() ; i++)
		{
			// sub j counter
			int sub_j = 0;
			// for each zero (direction spanned by the sub-hyper-cube)
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (c.c[j] == 0)
				{
					// take the combination from the sub-hyper-cube
					vc[i].c[j] = v[i].c[sub_j];
					sub_j++;
				}
				else
				{
					// take the combination from the hyper-cube position
					vc[i].c[j] = c.c[j];
				}
			}
		}

		// return the combinations
		return vc;
	}
};

#endif
