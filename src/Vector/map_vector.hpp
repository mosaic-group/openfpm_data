/*
 * map_vector.hpp
 *
 *  Created on: Aug 30, 2014
 *      Author: Pietro Incardona
 */

#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include <iostream>
#include <typeinfo>
#include "util/common.hpp"
#include "memory/PtrMemory.hpp"
#include "util/object_util.hpp"
#include "Grid/util.hpp"
#include "Vector/util.hpp"
#include "Vector/map_vector_grow_p.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "util/util_debug.hpp"
#include "Pack_stat.hpp"
#include "Grid/map_grid.hpp"
#include "memory/HeapMemory.hpp"
#include "vect_isel.hpp"
#include "util/object_s_di.hpp"
#include "util.hpp"
#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace openfpm
{
	/*! \brief Vector iterator
	 *
	 */

	class vector_key_iterator
	{
		//! Linearized end element
		size_t end;

	protected:

		//! Actual key
		size_t gk;

	public:

		/*! \brief Constructor require the size of the vector
		 *
		 * \param end size of the vector
		 */
		vector_key_iterator(size_t end, size_t start = 0)
		: end(end),gk(start)
		{}


		/*! \brief Get the next element
		 *
		 * Get the next element
		 *
		 * \return the next grid_key
		 *
		 */
		vector_key_iterator operator++()
		{
			//! increment the first index

			gk++;

			return *this;
		}

		/*! \brief Set the dimension
		 *
		 * \param d is the dimension (IGNORED is by default 0)
		 * \param sz set the counter to sz
		 *
		 */
		void set(int d, size_t sz)
		{
			// set the counter dim to sz

			gk = sz;
		}

		/*! \brief Check if there is the next element
		 *
		 * Check if there is the next element
		 *
		 * \return true if there is the next, false otherwise
		 *
		 */
		bool isNext()
		{
			if (gk < end)
			{
				//! we did not reach the end of the grid

				return true;
			}

			//! we reach the end of the grid
			return false;
		}

		/*! \brief Get the actual key
		 *
		 * Get the actual key
		 *
		 * \return the actual key
		 *
		 */
		size_t get()
		{
			return gk;
		}
	};
	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * Stub object look at the various implementations
	 *
	 * \snippet vector_test_util.hpp Create add and access
	 *
	 * \param T type of structure the vector has to store
	 * \param Memory allocator to use
	 * \param grow_p grow policy, how this vector should grow
	 *
	 * \see vector<T,HeapMemory,grow_policy_double,STD_VECTOR>
	 * \see vector<T,Memory,grow_p,OPENFPM_NATIVE>
	 * \see vector<T,Memory,grow_p,OPENFPM_NATIVE>
	 *
	 */
	template<typename T, typename Memory, typename grow_p, unsigned int impl>
	class vector
	{
	};

	#include "map_vector_std.hpp"

	/*! \brief Implementation of 1-D std::vector like structure
	 *
	 * The layout is memory_traits_lin
	 *
	 * ### Add and access elements
	 * \snippet vector_test_util.hpp Create add and access
	 *
	 * \tparam T type of object the vector store
	 * \tparam base memory layout to use
	 * \tparam Memory allocator to use
	 * \tparam grow_p grow policy for vector in case of reallocation
	 *
	 * OPENFPM_NATIVE implementation
	 *
	 */
	template<typename T,typename Memory, typename grow_p>
	class vector<T,Memory,grow_p,OPENFPM_NATIVE>
	{
		//! This structure use this layout
		typedef typename grid_cpu<1,T,Memory>::memory_lin memory_lin;

		//! Actual size of the vector, warning: it is not the space allocated in grid
		//! grid size increase by a fixed amount every time we need a vector bigger than
		//! the actually allocated space
		size_t v_size;

		//! 1-D static grid
		grid_cpu<1,T,Memory> base;

		/*! \brief Create a 1D vector that contain the vector size
		 *
		 * Used to construct the underline 1D-grid
		 *
		 * \param sz size of the vector
		 *
		 * \return a vector with one element
		 *
		 */
		std::vector<size_t> getV(size_t sz)
		{
			std::vector<size_t> tmp;

			tmp.push_back(sz);

			return tmp;
		}

		/*! \brief If the argument is zero return 1 otherwise return the argument
		 *
		 * \param sz output
		 * \param arg argument
		 *
		 */
		void non_zero_one(size_t sz[1], size_t arg)
		{
			if (arg == 0)
			{sz[0] = 1;}
			else
			{sz[0] = arg;}
		}

	public:

		//! it define that it is a vector
		typedef int yes_i_am_vector;

		//! Type of the encapsulation memory parameter
		typedef typename grid_cpu<1,T>::memory_conf memory_conf;

		//! iterator for the vector
		typedef vector_key_iterator iterator_key;

		//! Object container for T, it is the return type of get_o it return a object type trough
		// you can access all the properties of T
		typedef typename grid_cpu<1,T>::container container;

		//! Type of the value the vector is storing
		typedef T value_type;

		/*! \brief pack a vector selecting the properties to pack
		 *
		 * \param mem preallocated memory where to pack the vector
		 * \param obj object to pack
		 * \param sts pack-stat info
		 *
		 */

		struct testtest
		{
			static bool noPointers()	{return true;}
		};


		template<int ... prp> void pack(ExtPreAlloc<Memory> & mem, openfpm::vector<T> & obj, Pack_stat & sts)
		{
	#ifdef DEBUG
			if (mem.ref() == 0)
				std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
	#endif

			// if no properties should be packed return
			if (sizeof...(prp) == 0)
				return;

			if (has_Pack<T>::type::value == true)
				for (int i = 0; i < obj.size(); i++) {
					T var = obj.get(i);
					var.pack<prp...>(mem, sts);
				}


			// Sending property object
			typedef openfpm::vector<T> vctr;
			typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;

			typedef openfpm::vector<prp_object,ExtPreAlloc<Memory>,openfpm::grow_policy_identity> dtype;

			// Create an object over the preallocated memory (No allocation is produced)
			dtype dest;
			dest.setMemory(mem);
			dest.resize(obj.size());
			auto obj_it = obj.getIterator();

			while (obj_it.isNext())
			{
				// copy all the object in the send buffer
				typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_src;
				// destination object type
				typedef encapc<1,prp_object,typename dtype::memory_conf > encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src,encap_dst,ENCAP,prp...>(obj.get(obj_it.get()),dest.get(obj_it.get()));

				++obj_it;
			}

			// Update statistic
			sts.incReq();

		}

		/*! \brief Insert an allocation request into the vector
		 *
		 * \param obj vector object to pack
		 * \param requests vector
		 *
		 */
		template<int ... prp> void packRequest(openfpm::vector<T> & obj, std::vector<size_t> & v)
		{
			typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

			// Calculate the required memory for packing
			size_t alloc_ele = openfpm::vector<prp_object>::calculateMem(obj.size(),0);

			v.push_back(alloc_ele);
		}

		/*! \brief unpack a vector
		 *
		 * \warning the properties should match the packed properties, and the obj must have the same size of the packed vector, consider to pack
		 *          this information if you cannot infer-it
		 *
		 * \param ext preallocated memory from where to unpack the vector
		 * \param obj object where to unpack
		 *
		 */
		template<unsigned int ... prp> void unpack(ExtPreAlloc<Memory> & mem, openfpm::vector<T> & obj, Unpack_stat & ps)
		{
			// if no properties should be unpacked return
			if (sizeof...(prp) == 0)
				return;

			size_t id = 0;

			// Sending property object
			typedef openfpm::vector<T> vctr;
			typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;
			typedef openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> stype;

			// Calculate the size to pack the object
			size_t size = stype::calculateMem(obj.size(),0);

			// Create a Pointer object over the preallocated memory (No allocation is produced)
			PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

			stype src;
			src.setMemory(ptr);
			src.resize(obj.size());
			auto obj_it = obj.getIterator();

			while (obj_it.isNext())
			{
				// copy all the object in the send buffer
				typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_dst;
				// destination object type
				typedef encapc<1,prp_object,typename stype::memory_conf > encap_src;

				// Copy only the selected properties
				object_s_di<encap_src,encap_dst,ENCAP,prp...>(src.get(id),obj.get(obj_it.get()));

				++id;
				++obj_it;
			}

			ps.addOffset(size);
		}


		/*! \brief Return the size of the vector
		 *
		 * \return the size
		 *
		 */
		size_t size() const
		{
			return v_size;
		}

		/*! \brief Reserve slots in the vector to avoid reallocation
		 *
		 * Reserve slots in the vector to avoid reallocation
		 *
		 * \param sp number of slot to reserve
		 *
		 */

		void reserve(size_t sp)
		{
			if (sp > base.size())
			{
				//! Resize the memory
				size_t sz[1] = {sp};
				base.template resize<Memory>(sz);
			}
		}

		/*! \brief Resize the vector
		 *
		 * Resize the vector and allocate n elements
		 *
		 * \param slot number of elements
		 *
		 */
		void resize(size_t slot)
		{
			// If we need more space than what we allocated, allocate new memory

			if (slot > base.size())
			{
				size_t gr = grow_p::grow(base.size(),slot);

				//! Resize the memory
				size_t sz[1] = {gr};
				base.resize(sz);
			}

			// update the vector size
			v_size = slot;
		}

		//! Access key for the vector
		typedef size_t access_key;

		/*! \brief It insert a new emtpy object on the vector, eventually it reallocate the grid
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add()
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the grid
		 *
		 * \param v element to add
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 */
		void add(const T & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! copy the element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief It insert a new object on the vector, eventually it reallocate the vector
		 *
		 * \param v object (encapsulated)
		 *
		 * \warning It is not thread safe should not be used in multi-thread environment
		 *          reallocation, work only on cpu
		 *
		 *
		 */
		void add(const typename grid_cpu<1,T>::container & v)
		{
			//! Check if we have enough space

			if (v_size >= base.size())
			{
				//! Resize the memory, double up the actual memory allocated for the vector
				size_t sz[1];
				non_zero_one(sz,2*base.size());
				base.resize(sz);
			}

			//! copy the added element
			base.set(v_size,v);

			//! increase the vector size
			v_size++;
		}

		/*! \brief It add the element of another vector to this vector
		 *
		 * \param v from where to take the vector
		 *
		 */
		template <typename M, typename gp> void add(const vector<T, M,gp,OPENFPM_NATIVE> & v)
		{
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
				add(v.get(i));
		}

		/*! \brief It add the element of a source vector to this vector
		 *
		 * The number of properties in the source vector must be smaller than the destination
		 * all the properties of S must be mapped so if S has 3 properties
		 * 3 numbers for args are required
		 *
		 * \tparam S Base object of the source vector
		 * \tparam M memory type of the source vector
		 * \tparam gp Grow policy of the source vector
		 * \tparam args one or more number that define which property to set-up
		 *
		 * \param v source vector
		 *
		 */
		template <typename S, typename M, typename gp, unsigned int ...args> void add_prp(const vector<S, M,gp,OPENFPM_NATIVE> & v)
		{
			//! Add the element of v
			for (size_t i = 0 ; i < v.size() ; i++)
			{
				// Add a new element
				add();

				// write the object in the last element
				object_s_di<decltype(v.get(i)),decltype(get(size()-1)),ENCAP,args...>(v.get(i),get(size()-1));
			}
		}

		/*! \brief Remove one entry from the vector
		 *
		 * \param key element to remove
		 *
		 */
		void remove(size_t key)
		{
			size_t d_k = key;
			size_t s_k = key + 1;

			// keys
			while (s_k < size())
			{
				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size--;
		}

		/*! \brief Remove several entries from the vector
		 *
		 * \param keys objects id to remove
		 * \param start key starting point
		 *
		 */
		void remove(openfpm::vector<size_t> keys, size_t start = 0)
		{
			// Nothing to remove return
			if (keys.size() <= start )
				return;

			size_t a_key = start;
			size_t d_k = keys.get(a_key);
			size_t s_k = keys.get(a_key) + 1;

			// keys
			while (s_k < size())
			{
				// s_k should always point to a key that is not going to be deleted
				while (a_key+1 < keys.size() && s_k == keys.get(a_key+1))
				{
					a_key++;
					s_k = keys.get(a_key) + 1;
				}

				// In case of overflow
				if (s_k >= size())
					break;

				set(d_k,get(s_k));
				d_k++;
				s_k++;
			}

			// re-calculate the vector size

			v_size -= keys.size() - start;
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */

		template <unsigned int p>inline const typename type_cpu_prop<p,memory_lin>::type & get(size_t id) const
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " overflow" <<  "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline const typename grid_cpu<1,T>::container get(size_t id) const
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error: "  << __FILE__ << ":" << __LINE__ << " vector overflow" << "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \tparam p Property to get
		 * \param id Element to get
		 *
		 * \return the element value requested
		 *
		 */

		template <unsigned int p>inline typename type_cpu_prop<p,memory_lin>::type & get(size_t id)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " overflow" <<  "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.template get<p>(key);
		}

		/*! \brief Get an element of the vector
		 *
		 * Get an element of the vector
		 *
		 * \param id Element to get
		 *
		 * \return the element (encapsulated)
		 *
		 */

		inline typename grid_cpu<1,T>::container get(size_t id)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error "  << __FILE__ << ":" << __LINE__ << " vector overflow" << "\n";}
#endif
			grid_key_dx<1> key(id);

			return base.get_o(key);
		}

		/*! \brief Constructor from vector
		 *
		 * \param mv vector
		 *
		 */
		vector(vector<T, Memory,grow_p,OPENFPM_NATIVE> && mv)
		:v_size(mv.v_size),base(mv.base)
		{
		}

		/*! \brief It duplicate the vector
		 *
		 * \return a duplicated vector
		 *
		 */
		vector<T, Memory,grow_p,OPENFPM_NATIVE> duplicate() const
		{
			vector<T, Memory,grow_p,OPENFPM_NATIVE> dup;

			dup.v_size = v_size;
			dup.base.swap(base.duplicate());

			return dup;
		}

		/*! \brief Constructor from another temporal vector
		 *
		 * \param v the vector
		 *
		 */
		vector(const vector<T, Memory,grow_p,OPENFPM_NATIVE> && v)
		:v_size(0)
		{
			swap(v);
		}

		/*! \brief Constructor from another constant vector
		 *
		 * \param v the vector
		 *
		 */
		vector(const vector<T, Memory,grow_p,OPENFPM_NATIVE> & v)
		:v_size(0)
		{
			swap(v.duplicate());
		}

		//! Constructor, vector of size 0
		vector():v_size(0),base(getV(0))
		{
			base.setMemory();
		}

		//! Constructor, vector of size sz
		vector(size_t sz):v_size(sz),base(getV(sz))
		{
			base.setMemory();
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id element
		 * \param obj object (encapsulated)
		 *
		 */
		void set(size_t id, const typename grid_cpu<1,T>::container & obj)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the object id to obj
		 *
		 * \param id
		 * \param obj
		 *
		 */
		void set(size_t id, T & obj)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			//! copy the element
			base.set(id,obj);
		}

		/*! \brief Set the element of the vector v from another element of another vector
		 *
		 * \param id element id
		 * \param v vector source
		 * \param src source element
		 *
		 */
		void set(size_t id, vector<T,Memory,grow_p,OPENFPM_NATIVE> & v, size_t src)
		{
#ifdef DEBUG
			if (id >= v_size)
			{std::cerr << "Error " << __FILE__ << "  " << __LINE__ << " id overflow your vector" << "\n";}
#endif
			base.set(id,v.base,src);
		}

		/*! \brief Assignment operator
		 *
		 * move semantic movement operator=
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T, Memory,grow_p,OPENFPM_NATIVE> operator=(vector<T, Memory,grow_p,OPENFPM_NATIVE> && mv)
		{
			v_size = mv.v_size;
			base.swap(mv.base);

			return *this;
		}

		/*! \brief Assignment operator
		 *
		 * it copy
		 *
		 * \param mv vector
		 *
		 * \return itself
		 *
		 */
		vector<T, Memory,grow_p,OPENFPM_NATIVE> operator=(const vector<T, Memory,grow_p,OPENFPM_NATIVE> & mv)
		{
			v_size = mv.v_size;
			size_t rsz[1] = {v_size};
			base.resize(rsz);

			// copy the object
			for (size_t i = 0 ; i < v_size ; i++ )
			{
				grid_key_dx<1> key(i);
				base.set(key,mv.base,key);
			}

			return *this;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,Memory,grow_p,OPENFPM_NATIVE> & v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;
		}

		/*! \brief Swap the memory with another vector
		 *
		 * \param v vector
		 *
		 */
		void swap(openfpm::vector<T,Memory,grow_p,OPENFPM_NATIVE> && v)
		{
			size_t sz_sp = v_size;

			// swap the v_size
			v_size = v.v_size;

			base.swap(v.base);
			v.v_size = sz_sp;
		}

		/*! \brief Get iterator over the particles from a particular index
		 *
		 * \return an iterator to iterate from a particular index
		 *
		 */
		vector_key_iterator getIteratorFrom(size_t mark) const
		{
			return vector_key_iterator(v_size,mark);
		}

		/*! \brief Get the vector elements iterator
		 *
		 *
		 * \return an iterator to iterate through all the elements of the vector
		 *
		 */

		vector_key_iterator getIterator() const
		{
			return vector_key_iterator(v_size);
		}

		/*! \brief Return the size of the message needed to pack this object
		 *
		 * \return The size
		 *
		 */

		size_t packObjectSize()
		{
			return base.packObjectSize();
		}

		/*! \brief Pack the object into the given pointer
		 *
		 * \param mem pointer
		 *
		 * \return the size of the packed message
		 *
		 */
		size_t packObject(void * mem)
		{
			return base.packObject(mem);
		}

		/*! \brief Calculate the memory size required to allocate n elements
		 *
		 * Calculate the total size required to store n-elements in a vector
		 *
		 * \param n number of elements
		 * \param e unused
		 *
		 * \return the size of the allocation number e
		 *
		 */
		inline static size_t calculateMem(size_t n, size_t e)
		{
			if (n == 0)
				return 0;
			else
				return grow_p::grow(0,n) * sizeof(T);
		}

		/*! \brief How many allocation are required to create n-elements
		 *
		 * \param n number of elements
		 *
		 * \return the number of allocations
		 *
		 */
		inline static size_t calculateNMem(size_t n)
		{
			return 1;
		}

		/*! \brief Set the memory of the base structure using an object
		 *
		 * \param mem Memory object to use for allocation
		 *
		 */
		void setMemory(Memory & mem)
		{
			base.setMemory(mem);
		}

		/*! \brief Return the pointer that store the data
		 *
		 * \return the pointer that store the data
		 *
		 */
		void * getPointer()
		{
			return base.getPointer();
		}

		/*! \brief This class has pointer inside
		 *
		 * \return false
		 *
		 */
		static bool noPointers()
		{
			return false;
		}
	};

	template <typename T> using vector_std = vector<T, HeapMemory, openfpm::grow_policy_double, STD_VECTOR>;

}

#endif
