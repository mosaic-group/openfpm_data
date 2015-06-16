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

// Shift transformation
template<unsigned int dim, typename T>
class shift
{
	Point<dim,T> sh;

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
	}

	/*! \brief Shift the point transformation
	 *
	 * \param s source point
	 * \param i coordinate
	 *
	 * \return the transformed coordinate
	 *
	 */
	inline T transform(const T(&s)[dim], const size_t i)
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
	inline T transform(const Point<dim,T> & s, const size_t i)
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
	template<typename Mem> inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const size_t i)
	{
		return s.get(i) - sh.get(i);
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 *
	 */
	inline void setTransform(Matrix<dim,T> & mat, Point<dim,T> & orig)
	{
		for (size_t i = 0 ; i < dim ; i++)
			sh.get(i) = orig.get(i);
	}
};

// No transformation
template<unsigned int dim, typename T>
class no_transform
{
public:

	/*! \brief Constructor
	 *
	 * \param t Matrix transformation
	 * \param s shift
	 *
	 */
	no_transform(const Matrix<dim,T> & t, const Point<dim,T> & s)
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
	inline T transform(const T(&s)[dim], const size_t i)
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
	inline T transform(const Point<dim,T> & s, const size_t i)
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
	template<typename Mem> inline T transform(const encapc<1,Point<dim,T>,Mem> & s, const size_t i)
	{
		return s.template get<Point<dim,T>::x>()[i];
	}

	/*! \brief Set the transformation Matrix and shift
	 *
	 * \param mat Matrix transformation
	 * \param orig origin point
	 * \param crd coordinate
	 *
	 */
	inline void setTransform(Matrix<dim,T> & mat, Point<dim,T> & orig, size_t crd)
	{

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
 * \snippet CellList_test.hpp Cell decomposer use without shift
 * ### Cell decomposer with padding
 * \snippet CellList_test.hpp Test Cell decomposer with padding
 * ### Cell decomposer with shift
 * \snippet CellList_test.hpp Test Cell decomposer with shift
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

protected:

	// Total number of cell
	size_t tot_n_cell;

	// Domain of the cell list
	SpaceBox<dim,T> box;

	// Unit box of the Cell list
	SpaceBox<dim,T> box_unit;

	// Grid structure of the Cell list
	grid_sm<dim,void> gr_cell;

	// cell padding on each dimension
	size_t off[dim];

	/*! \brief Initialize all the structures
	 *
	 */
	void Initialize(const size_t pad, const size_t (& div)[dim])
	{
		// Space
		typedef Point<dim,T> p;

		// created a padded div
		size_t div_p[dim];

		for (size_t i = 0 ; i < dim ; i++)
			div_p[i] = div[i] + 2*pad;

		gr_cell.setDimensions(div_p);

		tot_n_cell = 1;

		// Total number of cells and calculate the unit cell size

		for (size_t i = 0 ; i < dim ; i++)
		{
			tot_n_cell *= gr_cell.size(i);

			// Cell are padded by
			box_unit.setHigh(i,box.getHigh(i) / (gr_cell.size(i)- 2*pad) );
		}

		for (size_t i = 0; i < dim ; i++)
			off[i] = pad;

		// Initialize p_middle

		p_middle = box_unit.getP2();
		p_middle = p_middle / 2;
	}

public:

	/*! \brief Return the underlying grid information of the cell list
	 *
	 * \return the grid infos
	 *
	 */
	grid_sm<dim,void> & getGrid()
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
	inline grid_key_dx<dim> getCellGrid(const T (& pos)[dim])
	{
#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif



		grid_key_dx<dim> key;
		key.set_d(0,t.tranform(pos[0]) / box_unit.getHigh(0) + off[0]);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef DEBUG
			if ((size_t)(t.transform(pos[s]) / box_unit.getHigh(s)) + off[s] < 0)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space";
#endif
			key.set_d(s,(size_t)(t.transform(pos[s]) / box_unit.getHigh(s)) + off[s]);
		}

		return key;
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
	grid_key_dx<dim> getCellGrid(const Point<dim,T> pos)
	{
#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		grid_key_dx<dim> key;
		key.set_d(0,t.transform(pos,0) / box_unit.getHigh(0) + off[0]);

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef DEBUG
			if ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s] < 0)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space";
#endif
			key.set_d(s,(size_t)(t.transform(pos,s) / box_unit.getHigh(s) + off[s]));
		}

		return key;
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
	size_t getCell(const T (& pos)[dim])
	{
#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		size_t cell_id = t.transform(pos,0) / box_unit.getHigh(0) + off[0];

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef DEBUG
			if (((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s]) < 0)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space";
#endif
			cell_id += gr_cell.size(s) * ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s]);
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
	size_t getCell(const Point<dim,T> & pos)
	{
#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif

		size_t cell_id = (size_t)(t.transform(pos,0) / box_unit.getHigh(0)) + off[0];

		for (size_t s = 1 ; s < dim ; s++)
		{
#ifdef DEBUG
			if (((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s]) < 0)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " point is not inside the cell space";
#endif
			cell_id += gr_cell.size_s(s-1) * ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s]);
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
	template<typename Mem> size_t getCell(const encapc<1,Point<dim,T>,Mem> & pos)
	{

#ifdef DEBUG
		if (tot_n_cell == 0)
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " using an uninitialized CellDecomposer";
#endif
		typedef Point<dim,T> p;

		size_t cell_id = (size_t)(t.transform(pos,0) / box_unit.getHigh(0)) + off[0];

		for (size_t s = 1 ; s < dim ; s++)
		{
			cell_id += gr_cell.size_s(s-1) * ((size_t)(t.transform(pos,s) / box_unit.getHigh(s)) + off[s]);
		}

		return cell_id;
	}

	/*! \brief Return the box smallest box containing the grid points
	 *
	 * Suppose a grid 5x5 defined on a Box<2,float> box({0.0,0.0},{1.0,1.0})
	 * and a Box defined Box<2,float>({0.4,0.4},{0.8,0.8}), it will return
	 * a Box<2,size_t> (2,2) and (4,4). A visualization of it is shown in the
	 * picture below.
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
	 * The big box is given by the object itself
	 *
	 * \return the box containing the grid points
	 *
	 */
	Box<dim,size_t> getGridPoints(Box<dim,T> & s_box)
	{
		// Box with inside grid
		Box<dim,size_t> bx;

		// Point p2
		Point<dim,T> p2 = s_box.getP2();
		p2 = p2 - p_middle;

		// Point p1
		Point<dim,T> p1 = s_box.getP1();
		p1 = p1 - p_middle;

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
	void setDimensions(SpaceBox<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad,div);
	}

	/*! \brief Set the domain to decompose
	 *
	 * \param box Domain to decompose
	 * \param div array with the number of cells on each dimensions
	 * \param pad padding cell
	 *
	 */
	void setDimensions(Box<dim,T> & box, const size_t (&div)[dim], const size_t pad)
	{
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad,div);
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
	void setDimensions(SpaceBox<dim,T> & box, const size_t (&div)[dim], Matrix<dim,T> & mat, Point<dim,T> & orig, const size_t pad)
	{
		t.setTransform(mat,orig);
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad);
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
	void setDimensions(Box<dim,T> & box, const size_t (&div)[dim], Matrix<dim,T> & mat, Point<dim,T> & orig, const size_t pad)
	{
		t.setTransform(mat,orig);
		this->box = box;
		this->gr_cell.setDimensions(div);
		Initialize(pad);
	}

	/*! \brief Constructor
	 *
	 * \param box Space where is defined the cell list (it is assumed p1 = {0, .... 0})
	 * \param div Reference array to the number of divisions on each dimensions
	 * \param mat Transformation matrix, the point is transformed as p' = mat * p
	 * \param orig origin of the cell decomposition
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
	CellDecomposer_sm(SpaceBox<dim,T> & box, size_t (&div)[dim], Matrix<dim,T> & mat, Point<dim,T> & orig, const size_t pad)
	:t(Matrix<dim,T>::identity(),Point<dim,T>::zero()),box(box),gr_cell()
	{
		Initialize(pad);
	}

	/*! \brief Constructor
	 *
	 * \param box Space where is defined the cell list (it is assumed p1 = {0, .... 0})
	 * \param div Reference array to the number of divisions on each dimensions
	 * \param orig origin of the cell decomposition
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
	CellDecomposer_sm(SpaceBox<dim,T> & box, size_t (&div)[dim], Point<dim,T> & orig, const size_t pad)
	:t(Matrix<dim,T>::identity(),orig),box(box),gr_cell(div)
	{
		Initialize(pad,div);
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
	 * \endverbatim
	 *
	 * Cell with p are padding cell cell that are around but external the box, the cell number 9 that
	 * is at the origin of the box is identified with 9
	 *
	 */
	CellDecomposer_sm(SpaceBox<dim,T> & box, size_t (&div)[dim], const size_t pad)
	:t(Matrix<dim,T>::identity(),Point<dim,T>::zero_p()),box(box),gr_cell()
	{
		Initialize(pad,div);
	}


	//! default constructor
	CellDecomposer_sm()
	:t(Matrix<dim,T>::identity(),Point<dim,T>::zero_p()),tot_n_cell(0)
	{

	}
};


#endif /* CELLDECOMPOSER_HPP_ */
