/*! \brief This class read a Grid raw data format
 *
 * A grid row data format is very simple. The first n numbers indicate the
 * size in every dimension of the grid. The other is the data contained by the grid.
 * In particular if we are in 3D and we are saving a 45x50x30 grid after the 3 numbers
 * I am expecting  45x50x30 = 67500 objects of type T. There is no check the dimensionality
 * is correct, neither the type is correct
 *
 * \tparam dim dimensionality of the grid
 * \tparam T type of the grid
 *
 */
#include <iostream>
#include <Grid/map_grid.hpp>
#include <fstream>

#define FORTRAN_STYLE 1
#define STRUCT_OF_ARRAY 2

/*! \brief This is the scalar case
 *
 * \tparam T scalar type
 *
 */
template<unsigned int dim, typename Tg, typename Tr, unsigned int i>
struct meta_raw_read
{
	static inline void read(grid_cpu<dim,Tg> & gr,std::ifstream & raw)
	{
		auto it = gr.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			raw.read((char *)&gr.template get<i>(key),sizeof(Tr));

			++it;
		}
	}
};

/*! \brief This is the vector case
 *
 * \tparam T vector type
 *
 */
template<unsigned int dim, typename Tg, typename Tr ,unsigned int i, unsigned int nv>
struct meta_raw_read<dim,Tg,Tr[nv],i>
{
	static inline void read(grid_cpu<dim,Tg> & gr,std::ifstream & raw)
	{
		for (size_t k = 0 ; k < nv ; k++)
		{
			auto it = gr.getIterator();

			while (it.isNext())
			{
				auto key = it.get();

				raw.read((char *)&gr.template get<i>(key)[k],sizeof(Tr));

				++it;
			}
		}
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to read each property
 *
 * \tparam ele_g element that store the grid and its attributes
 * \param St type of space where the grid live
 *
 */
template<unsigned int dim, typename Tg>
struct raw_read
{
	//! Grid out
	grid_cpu<dim,Tg> & gr;

	//! File stream
	std::ifstream & fl;

	/*! \brief constructor
	 *
	 * \param gr grid to fill
	 * \param fl file from where to read
	 *
	 */
	raw_read(grid_cpu<dim,Tg> & gr,std::ifstream & fl)
	:gr(gr),fl(fl)
	{};

	//! It read for each property
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::mpl::at<typename Tg::type,boost::mpl::int_<T::value>>::type Tr;

    	meta_raw_read<dim,Tg,Tr,T::value>::read(gr,fl);
    }
};

template <unsigned int dim, typename T, typename idx_type>
class GridRawReader
{
public:

	//! Constructor
	GridRawReader()	{};


	/*! \brief Read a raw grid
	 *
	 * \param file raw file to read
	 * \param gr grid to fill
	 * \param opt option (FORTRAN_STYLE)
	 * \param skip skip N byte
	 *
	 */
	bool read(std::string file, grid_cpu<dim,T> & gr, size_t opt = 0, size_t skip = 0)
	{
		idx_type tmp;
		std::ifstream raw;
		raw.open (file, std::ios::binary );

		if (raw.is_open() == false)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error in opening the file: " << file << std::endl;
			return false;
		}

		// get length of file:
		raw.seekg (0, std::ios::end);
		size_t length = raw.tellg();
		raw.seekg (skip, std::ios::beg);

		if (opt & FORTRAN_STYLE)
			raw.read((char *)&tmp,sizeof(idx_type));

		size_t sz[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{
			sz[i] = 0;
			raw.read((char *)&sz[i],sizeof(idx_type));
		}

		if (opt & FORTRAN_STYLE)
			raw.read((char *)&tmp,sizeof(idx_type));

		if (opt & FORTRAN_STYLE)
			raw.read((char *)&tmp,sizeof(idx_type));

		grid_sm<dim,void> gs(sz);

		size_t offset = 0;
		if (opt & FORTRAN_STYLE)
			offset += 2*sizeof(idx_type)*2;

		// Check the the file size make sense
		if (length - dim*sizeof(idx_type) - offset - skip != gs.size()*sizeof(T) )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " Error the size of the raw file does not match with the calculated one" << std::endl;
			return false;
		}

		gr.setMemory();

		// resize the grid
		gr.resize(sz);

		if (!(opt & STRUCT_OF_ARRAY))
		{
			// read the data
			raw.read((char *)gr.getPointer(),gr.size()*sizeof(T));
			raw.close();
		}
		else
		{
			// for each property
			raw_read<dim,T> rr(gr,raw);

			boost::mpl::for_each< boost::mpl::range_c<int,0, T::max_prop> >(rr);
		}

		return true;
	}
};
