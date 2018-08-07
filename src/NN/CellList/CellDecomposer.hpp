/*
 * CellDecomposer.hpp
 *
 *  Created on: Apr 1, 2015
 *      Author: Pietro Incardona
 */

#ifndef CELLDECOMPOSER_HPP_
#define CELLDECOMPOSER_HPP_

#include "Space/SpaceBox.hpp"
#include "Space/Matrix.hpp"
#include "util/copy_compare/meta_compare.hpp"
#include "Grid/grid_sm.hpp"

#define CELL_DECOMPOSER 8001lu




/*! This Class apply a shift transformation before converting to Cell-ID
 *
 *
 */
template<unsigned int dim, typename T>
class shift
{
	//! Shift point
	Point<dim,T> sh;

	//! Matrix transformation
	Matrix<dim,T> mat;

public:

	/*! \brief Constructor
	 *
	 * \param t Matrix transformation
	 * \param s shift
	 *
	 */
	shift(const Matrix<dim,T> & t, const Point<dim,T> & s)
	:sh(s)
	{
		mat.identity();
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	inline T transform(const T(&s)[dim], const size_t i) const
	{
		return s[i] - sh.get(i);
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	inline T transform(const Point<dim,T> & s, const size_t i) const
	{
		return s.get(i) - sh.get(i);
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	template<typename Mem> inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const size_t i) const
	{
		return s.template get<0>()[i] - sh.get(i);
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 *
	 */
	inline void setTransform(const Matrix<dim,T> & mat, const Point<dim,T> & orig)
	{
		for (size_t i = 0 ; i < dim ; i++)
			sh.get(i) = orig.get(i);
	}

	/*! \brief Get the shift vector
	 *
	 * \return the shift vector
	 *
	 */
	inline const Point<dim,T> & getOrig() const
	{
		return sh;
	}

	/*! \brief Get the transformation Matrix
	 *
	 * \return the transformation matrix
	 *
	 */
	inline const Matrix<dim,T> & getMat() const
	{
		return mat;
	}

	/*! \brief It return true if the shift match
	 *
	 * \param s shift to compare with
	 *
	 * \return true if it match
	 *
	 */
	inline bool operator==(const shift<dim,T> & s)
	{
		return sh == s.sh;
	}

	/*! \brief It return true if the shift is different
	 *
	 * \param s shift to compare with
	 *
	 * \return true if the shift is different
	 *
	 */
	inline bool operator!=(const shift<dim,T> & s)
	{
		return !this->operator==(s);
	}
};

/*! This Class apply a shift transformation before converting to Cell-ID
 *
 *
 */
template<unsigned int dim, typename T>
class shift_only
{
	//! Shift point
	Point<dim,T> sh;

public:

	/*! \brief Constructor
	 *
	 * \param t Matrix transformation
	 * \param s shift
	 *
	 */
	shift_only(const Matrix<dim,T> & t, const Point<dim,T> & s)
	:sh(s)
	{
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	__device__ __host__  inline T transform(const T * s, const int i) const
	{
		return s[i] - sh.get(i);
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	__device__ __host__  inline T transform(const T(&s)[dim], const int i) const
	{
		return s[i] - sh.get(i);
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	__device__ __host__  inline T transform(const Point<dim,T> & s, const int i) const
	{
		return s.get(i) - sh.get(i);
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	template<typename Mem> __device__ __host__   inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const int i) const
	{
		return s.template get<0>()[i] - sh.get(i);
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 *
	 */
	inline void setTransform(const Matrix<dim,T> & mat, const Point<dim,T> & orig)
	{
		for (size_t i = 0 ; i < dim ; i++)
			sh.get(i) = orig.get(i);
	}

	/*! \brief Get the shift vector
	 *
	 * \return the shift vector
	 *
	 */
	inline const Point<dim,T> & getOrig() const
	{
		return sh;
	}


	/*! \brief It return true if the shift match
	 *
	 * \param s shift to compare with
	 *
	 * \return true if it match
	 *
	 */
	inline bool operator==(const shift<dim,T> & s)
	{
		return sh == s.sh;
	}

	/*! \brief It return true if the shift is different
	 *
	 * \param s shift to compare with
	 *
	 * \return true if the shift is different
	 *
	 */
	inline bool operator!=(const shift<dim,T> & s)
	{
		return !this->operator==(s);
	}
};

//! No transformation
template<unsigned int dim, typename T>
class no_transform
{
	//! shift transform
	Point<dim,T> orig;

	//! Matrix transform
	Matrix<dim,T> mat;

public:

	/*! \brief Constructor
	 *
	 * \param t Matrix transformation
	 * \param s shift
	 *
	 */
	no_transform(const Matrix<dim,T> & t, const Point<dim,T> & s)
	{
		orig.zero();
		mat.identity();
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	inline T transform(const T(&s)[dim], const size_t i) const
	{
		return s[i];
	}

	/*! \brief No transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the source point coordinate
	 *
	 */
	inline T transform(const Point<dim,T> & s, const size_t i) const
	{
		return s.get(i);
	}

	/*! \brief No transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the point coordinate
	 *
	 */
	template<typename Mem> inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const size_t i) const
	{
		return s.template get<Point<dim,T>::x>()[i];
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 *
	 */
	inline void setTransform(const Matrix<dim,T> & mat, const Point<dim,T> & orig)
	{

	}

	/*! \brief It return always true true
	 *
	 * There is nothing to compare
	 *
	 * \param nt unused
	 *
	 * \return true
	 *
	 */
	inline bool operator==(const no_transform<dim,T> & nt)
	{
		return true;
	}

	/*! \brief It return always false
	 *
	 * There is nothing to compare they cannot be differents
	 *
	 * \param nt unused
	 *
	 * \return false
	 *
	 */
	inline bool operator!=(const no_transform<dim,T> & nt)
	{
		return false;
	}

	/*! \brief Get the origin
	 *
	 * \return the shift vector
	 *
	 */
	inline const Point<dim,T> & getOrig() const
	{
		return orig;
	}

	/*! \brief Get the transformation Matrix
	 *
	 * \return get the return the matrix
	 *
	 */
	inline const Matrix<dim,T> & getMat() const
	{
		return mat;
	}
};

//! No transformation
template<unsigned int dim, typename T>
class no_transform_only
{

public:

	/*! \brief Constructor
	 *
	 * \param t Matrix transformation
	 * \param s shift
	 *
	 */
	no_transform_only(const Matrix<dim,T> & t, const Point<dim,T> & s)
	{
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	__device__ __host__ inline T transform(const T * s, const int i) const
	{
		return s[i];
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	__device__ __host__ inline T transform(const T(&s)[dim], const int i) const
	{
		return s[i];
	}

	/*! \brief No transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the source point coordinate
	 *
	 */
	__device__ __host__  inline T transform(const Point<dim,T> & s, const int i) const
	{
		return s.get(i);
	}

	/*! \brief No transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the point coordinate
	 *
	 */
	template<typename Mem> __device__ __host__ inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const int i) const
	{
		return s.template get<Point<dim,T>::x>()[i];
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 *
	 */
	inline void setTransform(const Matrix<dim,T> & mat, const Point<dim,T> & orig)
	{

	}

	/*! \brief It return always true true
	 *
	 * There is nothing to compare
	 *
	 * \param nt unused
	 *
	 * \return true
	 *
	 */
	inline bool operator==(const no_transform<dim,T> & nt)
	{
		return true;
	}

	/*! \brief It return always false
	 *
	 * There is nothing to compare they cannot be differents
	 *
	 * \param nt unused
	 *
	 * \return false
	 *
	 */
	inline bool operator!=(const no_transform<dim,T> & nt)
	{
		return false;
	}
};

/*! \brief Decompose a space into cells
 *
 * It is a convenient class for cell decomposition of an N dimensional space into cells
 *  and index linearization with getCell. Transform parameter is used to shift the space
 *  from the origin where the cell list is defined. It basically apply a transformation to the
 *  point every time we call getCell of getCellGrid. The shift is given by the starting point (p1)
 *  of the box that define where the Cell decomposition live
 *
 * \tparam dim dimension of the space divided by the cell
 * \tparam T information the cell contain
 * \tparam transform type of transformation (no_transformation shift ... ), carefully if p1 it different from {0, ... 0}
 *         shift class must be used
 *
 * Example of a 2D space decompised into Cells, 6x6 structure with padding 1 without shift, cell indicated with p are padding cell
 * the origin of the cell or point (0,0) is marked with cell number 9
 *
 * \verbatim
 * +-----------------------+
 * |p |p |p |p |p |p |p |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |  |  |  |  |  |  |p |
 * +-----------------------+
 * |p |9 |  |  |  |  |  |p |
 * +-----------------------+
 * |p |p |p |p |p |p |p |p |
 * +-----------------------+
 * \endverbatim
 *
 * ### Cell decomposer without shift
 * \snippet CellDecomposer_unit_tests.hpp Cell decomposer use without shift
 * ### Cell decomposer with padding
 * \snippet CellDecomposer_unit_tests.hpp Test Cell decomposer with padding
 * ### Cell decomposer with shift
 * \snippet CellDecompose_unit_tests.hpp Test Cell decomposer with shift
 *
 */
template<unsigned int dim,typename T, typename transform = no_transform<dim,T>>
class CellDecomposer_sm
{
	// Point in the middle of a cell
	//
	// \verbatim
	//
	//         C (0.1,0.1)
	// +-----+
	// |     |
	// |  .P |
	// |     |
	// +-----+
    //
	// \endverbatim
	//
	//    C is the cell and P is the point inside the middle of the cell
	//    for example if the cell is (0.1,0.1) P is (0.05,0.05)
	//
	//
	Point<dim,T> p_middle;

	// Point transformation before get the Cell object (useful for example to shift the cell list)
	transform t;

	/*! \brief Convert the coordinates into id
	 *
	 * \param x coordinate
	 * \param s dimension
	 *
	 */
	inline size_t ConvertToID(const T (&x)[dim] ,size_t s) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}

	/*! \brief Convert the coordinates into id
	 *
	 * \param x point
	 * \param s dimension
	 *
	 */
	inline size_t ConvertToID(const Point<dim,T> & x ,size_t s, size_t sc = 0) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}

	/*! \brief Convert the coordinates into id without apply shift
	 *
	 * \param x coordinate
	 * \param s dimension
	 *
	 */
	inline size_t ConvertToID_ns(const T (&x)[dim] ,size_t s) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1):id;
		return id;
	}

	/*! \brief Convert the coordinates into id without apply shift
	 *
	 * \param x point
	 * \param s dimension
	 *
	 */
	inline size_t ConvertToID_ns(const Point<dim,T> & x ,size_t s, size_t sc = 0) const
	{
		size_t id = openfpm::math::size_t_floor(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1):id;
		return id;
	}

	/*! \brief Convert the coordinates into id
	 *
	 * \param x point
	 * \param s dimension
	 *
	 */
	template <typename Mem> inline size_t ConvertToID_(const encapc<1,Point<dim,T>,Mem> & x ,size_t s, size_t sc = 0) const
	{
		size_t id = (size_t)(t.transform(x,s) / box_unit.getHigh(s)) + off[s];
		id = (id >= gr_cell.size(s))?(gr_cell.size(s)-1-cell_shift.get(s)):id-cell_shift.get(s);
		return id;
	}

	/*! \Check if a particle is outside the domain of the cell-list
	 *
	 * \param pos position of the particle
	 * \param s coordinate to check
	 *
	 */
	template<typename Ele> inline void check_and_print_error(const Ele & pos ,size_t s) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (pos[0] < box.getLow(s) - off[s]*box_unit.getP2()[s] || pos[s] > box.getHigh(s) + off[s]*box_unit.getP2()[s])
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << Point<dim,T>(pos).toString() << " is not inside the cell space";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif
	}


	template<typename Ele> inline size_t getCellDom_impl(const Ele & pos) const
	{
		check_and_print_error(pos,0);

		size_t cell_id = ConvertToID_ns(pos,0);

		cell_id = (cell_id == gr_cell.size(0) - off[0])?gr_cell.size(0) - off[0] - 1:cell_id;
		cell_id = (cell_id == off[0]-1)?off[0]:cell_id;

		cell_id -= cell_shift.get(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			/* coverity[dead_error_begin] */
			check_and_print_error(pos,s);

			size_t cell_idt = ConvertToID_ns(pos,s);

			cell_idt = (cell_idt == gr_cell.size(s) - off[s])?gr_cell.size(s) - off[s] - 1:cell_idt;
			cell_idt = (cell_idt == off[s]-1)?off[s]:cell_idt;

			cell_idt -= cell_shift.get(s);

			cell_id += gr_cell2.size_s(s-1) * cell_idt;
		}

		return cell_id;
	}

	template<typename Ele> inline size_t getCellPad_impl(const Ele & pos) const
	{
		check_and_print_error(pos,0);

		size_t cell_id = ConvertToID_ns(pos,0);
		cell_id = (cell_id == off[0])?off[0]-1:cell_id;
		cell_id = (cell_id == gr_cell.size(0) - off[0] - 1)?gr_cell.size(0) - off[0]:cell_id;

		cell_id -= cell_shift.get(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
			check_and_print_error(pos,s);

			size_t cell_idt = ConvertToID_ns(pos,s);
			cell_idt = (cell_idt == off[s])?off[s]-1:cell_idt;
			cell_idt = (cell_idt == gr_cell.size(s) - off[s] - 1)?gr_cell.size(s) - off[s]:cell_idt;

			cell_idt -= cell_shift.get(s);

			cell_id += gr_cell2.size_s(s-1) * cell_idt;
		}

		return cell_id;
	}


protected:

	//! Total number of cell
	size_t tot_n_cell;

	//! Domain of the cell list
	SpaceBox<dim,T> box;

	//! Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	//! Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

	//! Grid structure of the internal Cell list
	grid_sm<dim,void> gr_cell2;

	//! Box in continuum for the gr_cell2
	Box<dim,T> box_gr_cell2;

	//! cell padding on each dimension
	size_t off[dim];

	//! cell_shift
	Point<dim,long int> cell_shift;

	// temporal buffer
	mutable size_t div_wp[dim];

	/*! \brief Initialize all the structures
	 *
	 */
	void Initialize(const size_t pad, const size_t (& div)[dim])
	{
#ifdef SE_CLASS1

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (div[i] == 0)
			{
				std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the number of cells on each dimension must be different from zero\n";
				ACTION_ON_ERROR(CELL_DECOMPOSER)
			}
		}

#endif

		// created a padded div
		size_t div_p[dim];

		for (size_t i = 0 ; i < dim ; i++)
			div_p[i] = div[i] + 2*pad;

		gr_cell.setDimensions(div_p);
		gr_cell2.setDimensions(div_p);

		tot_n_cell = 1;

		// Total number of cells and calculate the unit cell size

		for (size_t i = 0 ; i < dim ; i++)
		{
			tot_n_cell *= gr_cell.size(i);

			// Cell are padded by
			box_unit.setHigh(i,(box.getHigh(i) - box.getLow(i)) / (gr_cell.size(i)- 2*pad) );
		}

		for (size_t i = 0; i < dim ; i++)
			off[i] = pad;

		// Initialize p_middle

		p_middle = box_unit.getP2();
		p_middle = p_middle / 2;
	}

public:

	/*! \brief Return the transformation
	 *
	 * \return the transform
	 *
	 */
	const transform & getTransform()
	{
		return t;
	}

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	const grid_sm<dim,void> & getGrid() const
	{
#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		return gr_cell;
	}

	/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	inline grid_key_dx<dim> getCellGrid(const T (& pos)[dim]) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s] < 0)
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space\n";
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			key.set_d(s,ConvertToID(pos,s));

		}

		return key;
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param cell-id in grid units
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellLin(grid_key_dx<dim> && k) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (gr_cell2.size(0) < k.get(0) + off[0] - cell_shift.get(0))
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " cell coordinate 0 = " << k.get(0) + off[0] - cell_shift.get(0) << " bigger than cell space " << gr_cell2.size(0) << std::endl;
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		size_t cell_id = k.get(0) + off[0] - cell_shift.get(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if (gr_cell2.size(s) < k.get(s) + off[s] - cell_shift.get(s))
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " cell coordinate 0 = " << k.get(0) + off[0] - cell_shift.get(0) << " bigger than cell space " << gr_cell2.size(s) << std::endl;
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			cell_id += gr_cell2.size_s(s-1) * (k.get(s) + off[s] -cell_shift.get(s));
		}

		return cell_id;
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param cell-id in grid units
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellLin(const grid_key_dx<dim> & k) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (gr_cell2.size(0) < k.get(0) + off[0] - cell_shift.get(0))
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " cell coordinate 0 = " << k.get(0) + off[0] - cell_shift.get(0) << " bigger than cell space " << gr_cell2.size(0) << std::endl;
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		size_t cell_id = k.get(0) + off[0] -cell_shift.get(0);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if (gr_cell2.size(s) < k.get(s) + off[s] - cell_shift.get(s))
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " cell coordinate 0 = " << k.get(0) + off[0] - cell_shift.get(0) << " bigger than cell space " << gr_cell2.size(s) << std::endl;
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			cell_id += gr_cell2.size_s(s-1) * (k.get(s) + off[s] -cell_shift.get(s));
		}

		return cell_id;
	}

	/*! \brief Get the cell-ids
	 *
	 * Convert the point coordinates into the cell ids (Careful it include padding)
	 *
	 * \param pos Point position
	 *
	 * \return the cell-ids ad a grid_key_dx<dim>
	 *
	 */
	inline grid_key_dx<dim> getCellGrid(const Point<dim,T> & pos) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer" << std::endl;
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		grid_key_dx<dim> key;
		key.set_d(0,ConvertToID(pos,0));

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s] < 0)
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space" << std::endl;
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			/* coverity[dead_error_line] */
			key.set_d(s,ConvertToID(pos,s));
		}

		return key;
	}

	/*! \brief Get the cell-id enforcing that is NOT a cell from the padding
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \note this function is in general used to bypass round-off error
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellDom(const Point<dim,T> & pos) const
	{
		return getCellDom_impl<Point<dim,T>>(pos);
	}


	/*! \brief Get the cell-id enforcing that is NOT a cell from the padding
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \note this function is in general used to bypass round-off error
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellDom(const T (& pos)[dim]) const
	{
		return getCellDom_impl<T[dim]>(pos);
	}

	/*! \brief Get the cell-id enforcing that is from a padding cell
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \note this function is in general used to bypass round-off error
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellPad(const Point<dim,T> & pos) const
	{
		return getCellPad_impl<Point<dim,T>>(pos);
	}

	/*! \brief Get the cell-id enforcing that is from a padding cell
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \note this function is in general used to bypass round-off error
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCellPad(const T (& pos)[dim]) const
	{
		return getCellPad_impl<T[dim]>(pos);
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCell(const T (& pos)[dim]) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (pos[0] < box.getLow(0) - off[0]*box_unit.getP2()[0] || pos[0] > box.getHigh(0) + off[0]*box_unit.getP2()[0])
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << toPointString(pos) << " is not inside the cell space";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		size_t cell_id = ConvertToID(pos,0);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if (pos[s] < box.getLow(s) - off[s]*box_unit.getP2()[s] || pos[s] > box.getHigh(s) + off[s]*box_unit.getP2()[s])
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << toPointString(pos) << " is not inside the cell space";
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			cell_id += gr_cell2.size_s(s-1) * ConvertToID(pos,s);
		}

		return cell_id;
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	inline size_t getCell(const Point<dim,T> & pos) const
	{
#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (pos.get(0) < box.getLow(0) - off[0]*box_unit.getP2()[0] || pos.get(0) > box.getHigh(0) + off[0]*box_unit.getP2()[0])
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << pos.toPointString() << " is not inside the cell space";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		size_t cell_id = ConvertToID(pos,0);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1
			if (pos.get(s) < box.getLow(s) - off[s]*box_unit.getP2()[s] || pos.get(s) > box.getHigh(s) + off[s]*box_unit.getP2()[s])
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << pos.toPointString() << " is not inside the cell space";
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			/* coverity[dead_error_line] */
			cell_id += gr_cell2.size_s(s-1) * ConvertToID(pos,s);
		}

		return cell_id;
	}

	/*! \brief Get the cell-id
	 *
	 * Convert the point coordinates into the cell id
	 *
	 * \param pos Point position
	 *
	 * \return the cell-id
	 *
	 */
	template<typename Mem> inline size_t getCell(const encapc<1,Point<dim,T>,Mem> & pos) const
	{

#ifdef SE_CLASS1
		if (tot_n_cell == 0)
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}

		if (pos.template get<0>()[0] < box.getLow(0) - off[0]*box_unit.getP2()[0] || pos.template get<0>()[0] > box.getHigh(0) + off[0]*box_unit.getP2()[0])
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << toPointString(pos) << " is not inside the cell space " << box.toString() << std::endl;
			ACTION_ON_ERROR(CELL_DECOMPOSER);
		}
#endif

		size_t cell_id = ConvertToID_(pos,0);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef SE_CLASS1

			if (pos.template get<0>()[s] < box.getLow(s) - off[s]*box_unit.getP2()[s] || pos.template get<0>()[s] > box.getHigh(s) + off[s]*box_unit.getP2()[s])
			{
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point " << toPointString(pos) << " is not inside the cell space";
				ACTION_ON_ERROR(CELL_DECOMPOSER);
			}
#endif
			cell_id += gr_cell2.size_s(s-1) * ConvertToID_(pos,s);
		}

		return cell_id;
	}

	/*! \brief Return the smallest box containing the grid points
	 *
	 * Suppose a grid 5x5 defined on a Box<2,float> box({0.0,0.0},{1.0,1.0})
	 * and feeding to the function a Box<2,float>({0.4,0.4},{0.8,0.8}), it will return
	 * a Box<2,size_t> (2,2) and (3,3). A visualization it is shown in the
	 * picture below. (the grid points are centered on each cell)
	 *
	 * \verbatim
	 *
	   	+----------+
		|. . . . . |
		|          |
		|. . . . . |
		|   +---+  |
		|. .|. .|. |
		|   |   |  |
		|. .|. .|. |
		|   +---+  |
		|. . . . . |
		+----------+

		\endverbatim
	 *
	 *
	 * \div grid size on each dimension
	 * \box Small box in the picture
	 *
	 * The big box is defined by "this" box
	 *
	 * \return the box containing the grid points
	 *
	 */
	inline Box<dim,size_t> getGridPoints(const Box<dim,T> & s_box) const
	{
		// Box with inside grid
		Box<dim,size_t> bx;

		// Point p2
		Point<dim,T> p2 = s_box.getP2();
		p2 = p2 - p_middle;

		// Point p1
		Point<dim,T> p1 = s_box.getP1();
		p1 = p1 + p_middle;

		bx.setP2(getCellGrid(p2));
		bx.setP1(getCellGrid(p1));

		return bx;
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param pad padding cell
	 *
	 */
	inline void setDimensions(const Box<dim,T> & box, const size_t (&div)[dim], const size_t (&div2)[dim], const size_t pad, Point<dim,long int> cell_shift)
	{
		Matrix<dim,T> mat;
		mat.identity();
		t.setTransform(mat,box.getP1());
		this->box = box;

		Initialize(pad,div);

		// calculate the box for gr_cell2
/*		Box<dim,T> box_gr_cell2;

		for (size_t i = 0 ; i < dim ; i++)
		{
			box_gr_cell2.setLow(i,cell_shift.get(i) * box_unit.getHigh(i));
			box_gr_cell2.setHigh(i,(cell_shift.get(i) + div2[i]) * box_unit.getHigh(i));
		}*/

		size_t cells[dim];

		for (size_t i = 0 ; i < dim ; i++)
			cells[i] = div2[i] + 2*pad;

		gr_cell2.setDimensions(cells);

		for (size_t i = 0 ; i < dim ; i++)
			this->cell_shift.get(i) = cell_shift.get(i) - off[i];
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param pad padding cell
	 *
	 */
	inline void setDimensions(const Box<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	{
		Matrix<dim,T> mat;
		mat.identity();
		t.setTransform(mat,box.getP1());
		this->box = box;
		Initialize(pad,div);
		this->cell_shift = 0;
	}

	/*! \brief Set the cell decomposition parameters + the nested
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param div2 array with the number of cells on each dimension for the nested decomposer
	 * \param mat transformation matrix the cell space is transformed by p' = A * p
	 * \param orig origin of the cell decomposition
	 * \param pad padding cell
	 *
	 */
	inline void setDimensions(const Box<dim,T> & box, const size_t (&div)[dim], const size_t (&div2)[dim], const Matrix<dim,T> & mat, const size_t pad, Point<dim,long int> cell_shift)
	{
		t.setTransform(mat,box.getP1());
		this->box = box;

		Initialize(pad,div);

		// The nested cell is big div2 + 2*off

		size_t div_with_off[dim];

		for(size_t i = 0 ; i < dim ; i++)
			div_with_off[i] = div2[i] + 2*off[i];

		gr_cell2.setDimensions(div_with_off);

		for (size_t i = 0 ; i < dim ; i++)
			this->cell_shift.get(i) = cell_shift.get(i) - off[i];
	}

	/*! \brief Set the cell decomposition parameters
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param mat transformation matrix the cell space is transformed by p' = A * p
	 * \param orig origin of the cell decomposition
	 * \param pad padding cell
	 *
	 */
	inline void setDimensions(const Box<dim,T> & box, const size_t (&div)[dim], const Matrix<dim,T> & mat, const size_t pad)
	{
		t.setTransform(mat,box.getP1());
		this->box = box;
		Initialize(pad,div);
		this->cell_shift = 0;
	}


	/*! \brief Set a consistent cell decomposer
	 *
	 * When we want to create a new extended CellDecomposer consistent with the old one
	 * this function can be used to guarantee such costrain.
	 * With consistent we mean that for each cell of the old CellDecomposer exist
	 * a Cell in the new CellDecomposer that perfectly overlap
	 *
	 * \param cd OLD CellDecomposer
	 * \param cell_extension extension of the CellDecomposer in term of cells to add in each directions (like Ghost)
	 *
	 */
	inline void setDimensions(const CellDecomposer_sm<dim,T,transform> & cd, const Box<dim,size_t> & cell_extension)
	{
		this->cell_shift = 0;

		// Get the space transformation

		t.setTransform(cd.getMat(),cd.getOrig());

		// The domain is equivalent to the old one
		this->box = cd.box;

		// The padding must be calculated

		size_t pad = 0;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (pad < cell_extension.getLow(i))
				pad = cell_extension.getLow(i);
			else if (pad > cell_extension.getHigh(i))
				pad = cell_extension.getHigh(i);
		}

		// We have to give the old division

		size_t sz_div[dim];

		for (size_t i = 0 ; i < dim ; i++)
			sz_div[i] = cd.gr_cell.size(i) - 2*cd.off[i];

		Initialize(pad,sz_div);
	}

	/*! \brief Constructor
	 *
	 * \param box Space where is defined the cell list
	 * \param div Reference array to the number of divisions on each dimensions
	 * \param mat Transformation matrix, the point is transformed as p' = mat * p
	 * \param pad cell padding
	 *
	 *  Example for div = {6,6} and pad = 1
	 *
	 * \verbatim
	 * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |9 |  |  |  |  |  |p |
     * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
	 * \endverbatim
	 *
	 * Cell with p are padding cell cell that are around but external the box, the cell number 9 that
	 * is at the origin of the box is identified with 9
	 *
	 */
	CellDecomposer_sm(const SpaceBox<dim,T> & box, const size_t (&div)[dim], Matrix<dim,T> & mat, const size_t pad)
	:t(Matrix<dim,T>::identity(),box.getP1()),box(box),gr_cell()
	{
		Initialize(pad);
	}

	/*! \brief Constructor
	 *
	 * \param box Space where is defined the cell list (it is assumed p1 = {0, .... 0})
	 * \param div Reference array to the number of divisions on each dimensions
	 * \param pad cell padding
	 *
	 *  Example for div = {7,7} and pad = 1
	 *
	 * \verbatim
	 * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |  |  |  |  |  |  |p |
     * +-----------------------+
     * |p |9 |  |  |  |  |  |p |
     * +-----------------------+
     * |p |p |p |p |p |p |p |p |
     * +-----------------------+
	 *
	 * \endverbatim
	 *
	 * Cell with p are padding cell cell that are around but external the box, the cell number 9 that
	 * is at the origin of the box is identified with 9
	 *
	 */
	CellDecomposer_sm(const SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	:t(Matrix<dim,T>::identity(),box.getP1()),box(box),gr_cell(div)
	{
		Initialize(pad,div);
	}

	/*! \brief Constructor for a consistent CellDecomposer construction
	 *
	 * \param cd Old CellDecomposer
	 * \param ext Extension box
	 *
	 *
	 */
	CellDecomposer_sm(const CellDecomposer_sm<dim,T,transform> & cd, Box<dim,size_t> & ext)
	:t(Matrix<dim,T>::identity(),cd.getOrig())
	{
		setDimensions(cd,ext);
	}


	//! default constructor
	CellDecomposer_sm()
	:t(Matrix<dim,T>::identity(),Point<dim,T>::zero_p()),tot_n_cell(0)
	{

	}

	/*! \brief Return the box that represent the cell
	 *
	 * Can be considered the spacing between vertices of the cells
	 *
	 * \return the box
	 *
	 */
	inline const Box<dim,T> & getCellBox() const
	{
		return box_unit;
	}

	/*! \brief Get the transformation Matrix of the cell decomposer
	 *
	 * \return the transformation Matrix
	 *
	 */
	inline const Matrix<dim,T> & getMat() const
	{
		return t.getMat();
	}

	/*! \brief Get the origin of the cell decomposer
	 *
	 * \return the origin
	 *
	 */
	inline const Point<dim,T> & getOrig() const
	{
		return t.getOrig();
	}

	/*! \brief Convert a Box in the domain space into cell units (Negative contour included, positive contour excluded)
	 *
	 *  Given the following
	 *
	 * \verbatim
	 *
                      +-----+-----+-----+-----+-----+-----+ (1.0. 1.0) Domain box
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
Box "b"      <-----------------+  |     |   | |     |     |  Grid (7, 6)
(0.1 , 0.42)          |   | |     |     |   | |     |     |  Cells (6,5)
(0.64, 0.85)          +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
                      |   | |     |     |   | |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                    (0.0, 0.0)


    + = grid points

    \verbatim

    It return a Box with P1 = (0,2), P2 = (3,4)

	 *
	 * \param b Box in domain space
	 * \param bc boundary conditions
	 *
	 * \return Box in grid units, if P2 < P1 the box does not include any grid points
	 *
	 */
	inline Box<dim,long int> convertDomainSpaceIntoCellUnits(const Box<dim,T> & b_d, const size_t (& bc)[dim]) const
	{
		Box<dim,long int> g_box;
		Box<dim,T> b = b_d;
		b -= getOrig();

		// Convert b into grid units
		b /= getCellBox().getP2();

		// Considering that we are interested in a box decomposition of the space
		// where each box does not intersect any other boxes in the decomposition we include the negative
		// countour and exclude the positive one. So ceilP1 do the job for P1 while ceilP2 - 1
		// do the job for P2

		b.floorP1();
		b.ceilP2();

		g_box = b;

		// Translate the box by the offset

		for (size_t i = 0 ; i < dim ; i++)
		{
			g_box.setLow(i,g_box.getLow(i) + off[i]);
			g_box.setHigh(i,g_box.getHigh(i) + off[i]);
		}

		// on the other hand with non periodic boundary condition, the positive border of the
		// sub-domain at the edge of the domain must be included

		Point<dim,size_t> p_move;

		for (size_t i = 0 ; i < dim ; i++)
		{
			// we are at the positive border (We are assuming that there are not rounding error in the decomposition)
			if (b_d.getHigh(i) == box.getHigh(i))
			{
				if (bc[i] == NON_PERIODIC)
				{
					// Here the positive boundary is included
					g_box.setHigh(i,gr_cell.size(i) - off[i]);
				}
				else
				{
					// Carefull in periodic gr_cell is one bigger than the non-periodic
					// and the positive boundary is excluded
					g_box.setHigh(i,gr_cell.size(i)-1 - off[i]);
				}
			}


			if (b_d.getLow(i) == box.getHigh(i))
			{
				if (bc[i] == NON_PERIODIC)
				{
					// The instruction is the same but the meaning is different
					// for this reason there is anyway a branch
					// Here the border is not included
					g_box.setLow(i,gr_cell.size(i) - off[i]);
				}
				else
				{
					// Carefull in periodic gr_cell is one bigger than the non-periodic
					// Here the border is included
					g_box.setLow(i,gr_cell.size(i) - off[i]);
				}
			}

			/////////// Low boundary

			// we are at the positive border (We are assuming that there are not rounding error in the decomposition)
			/* coverity[copy_paste_error] */
			if (b_d.getHigh(i) == box.getLow(i))
					g_box.setHigh(i,off[i]);


			if (b_d.getLow(i) == box.getLow(i))
				g_box.setLow(i,off[i]);
		}

		return g_box;
	}


	/*! \brief Convert a Box in the domain space into grid units (Negative contour included, positive contour excluded)
	 *
	 *  Given the following
	 *
	 * \verbatim
	 *
                      +-----+-----+-----+-----+-----+-----+ (1.0. 1.0) Domain box
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
Box "b"      <-----------------+  |     |   | |     |     |  Grid (7, 6)
(0.1 , 0.42)          |   | |     |     |   | |     |     |
(0.64, 0.85)          +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
                      |   | |     |     |   | |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                    (0.0, 0.0)


    + = grid points

    \verbatim

    It return a Box with P1 = (1,3), P2 = (3,4)

	 *
	 * \param b Box in domain space
	 * \param bc boundary conditions
	 *
	 * \return Box in grid units, if P2 < P1 the box does not include any grid points
	 *
	 */
	inline Box<dim,long int> convertDomainSpaceIntoGridUnits(const Box<dim,T> & b_d, const size_t (& bc)[dim]) const
	{
		Box<dim,long int> g_box;
		Box<dim,T> b = b_d;
		b -= getOrig();

		// Convert b into grid units
		b /= getCellBox().getP2();

		// Considering that we are interested in a box decomposition of the space
		// where each box does not intersect any other boxes in the decomposition we include the negative
		// countour and exclude the positive one. So ceilP1 do the job for P1 while ceilP2 - 1
		// do the job for P2

		b.ceilP1();

		// (we do -1 later)
		b.ceilP2();
		for (size_t i = 0 ; i < dim ; i++)	{b.setHigh(i,b.getHigh(i)-1);}

		g_box = b;

		// on the other hand with non periodic boundary condition, the positive border of the
		// sub-domain at the edge of the domain must be included

		Point<dim,size_t> p_move;

		for (size_t i = 0 ; i < dim ; i++)
		{
			// we are at the positive border (We are assuming that there are not rounding error in the decomposition)
			if (b_d.getHigh(i) == box.getHigh(i))
			{
				if (bc[i] == NON_PERIODIC)
				{
					// Here the positive boundary is included
					g_box.setHigh(i,gr_cell.size(i) - off[i]);
				}
				else
				{
					// Carefull in periodic gr_cell is one bigger than the non-periodic
					// and the positive boundary is excluded
					g_box.setHigh(i,gr_cell.size(i)-1 - off[i]);
				}
			}


			if (b_d.getLow(i) == box.getHigh(i))
			{
				if (bc[i] == NON_PERIODIC)
				{
					// The instruction is the same but the meaning is different
					// for this reason there is anyway a branch
					// Here the border is not included
					g_box.setLow(i,gr_cell.size(i) - off[i]);
				}
				else
				{
					// Carefull in periodic gr_cell is one bigger than the non-periodic
					// Here the border is included
					g_box.setLow(i,gr_cell.size(i) - off[i]);
				}
			}

			/////////// Low boundary

			// we are at the positive border (We are assuming that there are not rounding error in the decomposition)
			/* coverity[copy_paste_error] */
			if (b_d.getHigh(i) == box.getLow(i))
					g_box.setHigh(i,off[i]);


			if (b_d.getLow(i) == box.getLow(i))
				g_box.setLow(i,off[i]);
		}

		return g_box;
	}

	inline Box<dim,long int> convertDomainSpaceIntoCellUnitsNear(const Box<dim,T> & b_d) const
	{
		Box<dim,long int> g_box;
		Box<dim,T> b = b_d;
		b -= getOrig();

		// Convert b into grid units
		b /= getCellBox().getP2();

		// Considering that we are interested in a box decomposition of the space
		// where each box does not intersect any other boxes in the decomposition we include the negative
		// countour and exclude the positive one. So ceilP1 do the job for P1 while ceilP2 - 1
		// do the job for P2

		for (size_t i = 0 ; i < dim ; i++)
		{
				b.setLow(i,b.getLow(i) + 0.125);
				b.setHigh(i,b.getHigh(i) + 0.125);
		}

		b.floorP1();
		b.floorP2();

		g_box = b;

		// Translate the box by the offset

		for (size_t i = 0 ; i < dim ; i++)
		{
				g_box.setLow(i,g_box.getLow(i) + off[i]);
				g_box.setHigh(i,g_box.getHigh(i) + off[i]);
		}


		return g_box;
	}

	/*! \brief Convert a Box in grid units into the domain space (Negative contour included, positive contour excluded)
	 *
	 *  Given the following
	 *
	 * \verbatim

                      +-----+-----+-----+-----+-----+-----+ (1.0. 1.0) Domain box
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
Box "b"      <-----------------+  |     |   | |     |     |  Grid (7, 6)
(0.1 , 0.42)          |   | |     |     |   | |     |     |
(0.64, 0.85)          +-----+-----+-----+-----+-----+-----+
                      |   | |     |     |   | |     |     |
                      |   | |     |     |   | |     |     |
                      |   +-----------------+ |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      |     |     |     |     |     |     |
                      +-----+-----+-----+-----+-----+-----+
                    (0.0, 0.0)


    + = grid points

    \verbatim

    It return a Box with P1 = (1,3), P2 = (3,4)

	 *
	 * \param b Box in grid units
	 *
	 * \return Box in domain space, if P2 < P1 the method return an invalid box
	 *
	 */
	inline Box<dim,T> convertCellUnitsIntoDomainSpace(const Box<dim,long int> & b_d) const
	{
		Box<dim,T> be;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if ((long int)gr_cell.size(i) - (long int)off[i] == b_d.getLow(i))
				be.setLow(i,box.getHigh(i));
			else if ((long int)off[i] == b_d.getLow(i))
				be.setLow(i,box.getLow(i));
			else
				be.setLow(i,(b_d.getLow(i) - off[i]) * box_unit.getP2()[i] + box.getLow(i));

			if ((long int)gr_cell.size(i) - (long int)off[i] == b_d.getHigh(i))
				be.setHigh(i,box.getHigh(i));
			else if ((long int)off[i] == b_d.getHigh(i))
				be.setHigh(i,box.getLow(i));
			else
				be.setHigh(i,(b_d.getHigh(i) - off[i]) * box_unit.getP2()[i] + box.getLow(i));
		}

		return be;
	}

	/*! \brief Return the number of divisions of the Cell Decomposer (including padding)
	 *
	 * \return the number of divisions
	 *
	 */
	const size_t (& getDiv() const)[dim]
	{
		return gr_cell.getSize();
	}


	/*! \brief Return the number of divisions of the Cell Decomposer (without padding)
	 *
	 * \return the number of divisions
	 *
	 */
	const size_t (& getDivWP() const)[dim]
	{
		for (size_t i = 0 ; i < dim ; i++)
		{div_wp[i] = gr_cell.getSize()[i] - 2*getPadding(i);}

		return div_wp;
	}

	/*! \brief Return the domain where the CellDecomposer is defined
	 *
	 * \return The domain of the remote machine
	 *
	 */
	const Box<dim,T> & getDomain() const
	{
		return box;
	}

	/*! \brief it swap the content of two Cell Decomposer
	 *
	 * \param cd CellDecomposer to swap with
	 *
	 */
	inline void swap(CellDecomposer_sm<dim,T,transform> & cd)
	{
		// swap all the members
		p_middle.swap(cd.p_middle);

		// Point transformation before get the Cell object (useful for example to shift the cell list)
		transform t_t = t;
		t = cd.t;
		cd.t = t_t;

		// Total number of cell
		size_t tot_n_cell_t = tot_n_cell;
		tot_n_cell = cd.tot_n_cell;
		cd.tot_n_cell = tot_n_cell_t;

		box.swap(cd.box);
		box_unit.swap(cd.box_unit);
		gr_cell.swap(cd.gr_cell);
		gr_cell2.swap(cd.gr_cell2);

		for (size_t i = 0 ; i < dim ; i++)
		{
			size_t off_t = off[i];
			off[i] = cd.off[i];
			cd.off[i] = off_t;

			size_t cs_t = cell_shift.get(i);
			cell_shift.get(i) = cd.cell_shift.get(i);
			cd.cell_shift.get(i) = cs_t;
		}
	}

	/*! \brief Check that the CellDecomposer is the same
	 *
	 * \param cd Cell decomposer
	 *
	 * \return true if the two CellDecomposer are the same
	 *
	 */
	inline bool operator==(const CellDecomposer_sm<dim,T,transform> & cd)
	{
		if (meta_compare<Point<dim,T>>::meta_compare_f(p_middle,cd.p_middle) == false)
			return false;

		if (t != cd.t)
			return false;

		if (tot_n_cell != cd.tot_n_cell)
			return false;

		if (box != cd.box)
			return false;

		if (box_unit != cd.box_unit)
			return false;

		if (gr_cell != cd.gr_cell)
			return false;

		if (gr_cell2 != cd.gr_cell2)
			return false;

		for (size_t i = 0 ; i < dim ; i++)
		{
			if (off[i] != cd.off[i])
				return false;

			if (cell_shift.get(i) != cd.cell_shift.get(i))
				return false;
		}

		return true;
	}

	/*! \brief Check that the CellDecomposer is the same
	 *
	 * \param cd Cell decomposer
	 *
	 * \return true if the two CellDecomposer is different
	 *
	 */
	inline bool operator!=(const CellDecomposer_sm<dim,T,transform> & cd)
	{
		return ! this->operator==(cd);
	}

	/*! \brief Return the number of padding cells of the Cell decomposer
	 *
	 * \param i dimension
	 *
	 * \return the number of padding cells
	 *
	 */
	size_t getPadding(size_t i) const
	{
		return off[i];
	}

	/*! \brief Return the number of padding cells of the Cell decomposer as an array
	 *
	 *
	 * \return the number of padding cells
	 *
	 */
	size_t (& getPadding())[dim]
	{
		return off;
	}

	/*! \brief Return the index of the first cell in the domain
	 *
	 * This function make sense if the CellDecomposer has been
	 * created with setDimensions with div and div2
	 *
	 * \return the first domain cell
	 *
	 */
	grid_key_dx<dim> getStartDomainCell() const
	{
		grid_key_dx<dim> key;

		for (size_t i = 0 ; i < dim ; i++)
		{
			key.set_d(i, cell_shift.get(i));
		}

		return key;
	}

	/*! \brief Return the index of the last cell in the domain
	 *
	 * This function make sense if the CellDecomposer has been
	 * created with setDimensions with div and div2
	 *
	 * \return the last domain cell
	 *
	 */
	grid_key_dx<dim> getStopDomainCell() const
	{
		grid_key_dx<dim> key;

		for (size_t i = 0 ; i < dim ; i++)
		{
			key.set_d(i,cell_shift.get(i) + gr_cell2.size(i) - 2*getPadding(i) - 1);
		}

		return key;
	}

	/*! \brief Return the cell shift
	 *
	 * \return the cell shifting
	 *
	 */
	grid_key_dx<dim> getShift() const
	{
		grid_key_dx<dim> k;

		for (size_t i = 0 ; i < dim ; i++)
			k.set_d(i,cell_shift.get(i));

		return k;
	}

	/*! \brief Return the grid structure of the internal cell list
	 *
	 * \return the grid information
	 *
	 */
	const grid_sm<dim,void> & getInternalGrid() const
	{
		return gr_cell2;
	}
};


#endif /* CELLDECOMPOSER_HPP_ */
