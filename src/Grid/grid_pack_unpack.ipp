/*!
 * This file contains the implemetation of packer and unpacker for grid
 * Created on: July 8, 2016
 *     Author: Yaroslav Zaluzhnyi
 */


//Functions to check if the packing object is complex
static bool pack()
{
	return false;
}

/*! \brief Reset the pack calculation
 * 
 * \note in this case does nothing
 *
 */
void packReset()
{}

/*! \brief Pack calculate
 * 
 * \note in this case does nothing
 *
 */
 template<unsigned int ... prp, typename context_type>
void packCalculate(size_t & req, const context_type & ctx)
{}

static bool packRequest()
{
	return false;
}

static bool packMem()
{
	return false;
}

//These structures do a packing of a simple (no "pack()" inside) object 

//With specified properties
template<bool sel, int ... prp>
struct pack_simple_cond
{
	static inline void pack(const grid_base_impl<dim,T,S,layout,layout_base> & obj, ExtPreAlloc<S> & mem, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif
		
		//Pack the size of a grid
		for (size_t i = 0; i < dim; i++)
		{
			Packer<size_t, S>::pack(mem,obj.getGrid().size(i),sts);
		}

		// Sending property object and vector
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>,typename layout_base<prp_object>::type,layout_base,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());

		auto it = obj.getIterator();

		// copy all the object in the send buffer
		typedef encapc<dims,value_type,layout > encap_src;
		// destination object type
		typedef encapc<1,prp_object,typename dtype::layout_type > encap_dst;

		pack_with_iterator<!is_contiguos<prp...>::type::value || has_pack_gen<prp_object>::value,
							dim,
							decltype(obj),
							   encap_src,
		 	 	 	 	 	   encap_dst,
							   typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,
							   decltype(it),
							   dtype,
							   prp...>::pack(obj,it,dest);

		// Update statistic
		sts.incReq();
	}
};

//Without specified properties
template<int ... prp>
struct pack_simple_cond<true, prp ...>
{
	static inline void pack(const grid_base_impl<dim,T,S,layout,layout_base> & obj, ExtPreAlloc<S> & mem, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif
		
		//Pack the size of a grid
		for (size_t i = 0; i < dim; i++)
		{
			Packer<size_t, S>::pack(mem,obj.getGrid().size(i),sts);
		}		

		// Sending property object
		typedef openfpm::vector<T,ExtPreAlloc<S>,typename layout_base<T>::type,layout_base,openfpm::grow_policy_identity> dtype;
		
		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());
	
		auto obj_it = obj.getIterator();
	
		size_t id = 0;
		
		while (obj_it.isNext())
		{
			// Copy
			dest.get(id).set(obj.get_o(obj_it.get()));
			
			++obj_it;
			++id;
		}
	
		// Update statistic
		sts.incReq();
	}
};

//These structures do an upacking of a simple object (no pack() inside)

//With specified properties
template<bool sel, int ... prp>
struct unpack_simple_cond
{
	static inline void unpack(grid_base_impl<dim,T,S,layout,layout_base> & obj , ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
		size_t dims[dim];
	
	    //Unpack a size of a source grid
		for (size_t i = 0; i < dim; i++)
		{
			size_t u2 = 0;
			Unpacker<size_t, S>::unpack(mem,u2,ps);
			dims[i] = u2;
		}
		
		//Resize a destination grid
		obj.resize(dims);
		obj.setMemory();
		
		// object that store the information in mem
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory,typename layout_base<prp_object>::type,layout_base,openfpm::grow_policy_identity> stype;

		// Calculate the size to pack the object
		size_t size = obj.packMem<prp...>(obj.size(),0);
		
		// Create an object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

		// Create an object over a pointer (No allocation is produced)
		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());

		auto it = obj.getIterator();

		// copy all the object in the send buffer
		typedef encapc<dim,grid_base_impl<dim,T,S,layout,layout_base>::value_type,layout > encap_dst;
		// destination object type
		typedef encapc<1,prp_object,typename memory_traits_lin<prp_object>::type > encap_src;


		unpack_with_iterator<dim,
							 decltype(obj),
							 encap_src,
							 encap_dst,
							 typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,
							 decltype(it),
							 stype,
							 prp...>::unpack(obj,it,src);

		ps.addOffset(size);
	}
};

//Without specified properties
template<int ... prp>
struct unpack_simple_cond<true, prp ...>
{
	static inline void unpack(grid_base_impl<dim,T,S,layout,layout_base> & obj , ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
		size_t dims[dim];
	
		//Unpack a size of a source grid
		for (size_t i = 0; i < dim; i++)
		{
			size_t u2 = 0;
			Unpacker<size_t, S>::unpack(mem,u2,ps);
			dims[i] = u2;
		}
		
		//Resize a destination grid
		obj.resize(dims);
		obj.setMemory();
		
		// Sending property object
		typedef openfpm::vector<T,PtrMemory,typename layout_base<T>::type,layout_base,openfpm::grow_policy_identity> stype;

		// Calculate the size to pack the object
		size_t size = obj.packMem<prp...>(obj.size(),0);
		
		// Create a Pointer object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));
		
		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());
		
		auto obj_it = obj.getIterator();
		
		size_t id = 0;
		
		while (obj_it.isNext())
		{
			// Copy
			obj.get_o(obj_it.get()).set(src.get(id));
		
			++id;
			++obj_it;
		}
		
		ps.addOffset(size);
	}
};
	
	/*! \brief Insert an allocation request into the vector
	 *
	 * \tparam prp list of properties
	 *
	 * \param v vector of allocation sequence
	 *
	 */
	template<int ... prp> inline void packRequest(size_t & req) const
	{
		//Pushback a "sizeof" of dimension sizes of the grid
		for (size_t i = 0; i < dim; i++)
		{
			req += sizeof(this->getGrid().size(i));
		}
		
		// If all of the aggregate properties do not have a "pack()" member	
		if (has_pack_agg<T,prp...>::result::value == false)
		{
			size_t alloc_ele = this->packMem<prp...>(this->size(),0);
			req += alloc_ele;
		}
		//If at least one property has "pack()"
		else
		{
			auto key_it = this->getIterator();

			while (key_it.isNext())
			{
				auto k = key_it.get();
				//Call a pack request
				call_aggregatePackRequest<decltype(this->get_o(k)),S,prp ... >::call_packRequest(this->get_o(k),req);
				
				++key_it;
			}
		}
	}

	/*! \brief pack a grid selecting the properties to pack
	 *
	 * \param mem preallocated memory where to pack the grid
	 * \param sts pack-stat info
	 *
	 */
	template<int ... prp> inline void pack(ExtPreAlloc<S> & mem, Pack_stat & sts) const
	{	
		//If all of the aggregate properties are simple (don't have "pack()" member)
		if (has_pack_agg<T,prp...>::result::value == false)
		{
			//Call a packer
			pack_simple_cond<sizeof...(prp) == 0,prp...>::pack(*this,mem,sts);
		}
		//If at least one property has a "pack()" member
		else
		{
			//Pack the size of a grid
			for (size_t i = 0; i < dim; i++)
			{
				Packer<size_t, S>::pack(mem,this->getGrid().size(i),sts);
			}
			
			auto key_it = this->getIterator();

			while (key_it.isNext())
			{
				auto k = key_it.get();
				//Call a pack request
				call_aggregatePack<decltype(this->get_o(k)),S,prp ... >::call_pack(this->get_o(k),mem,sts);
				
				++key_it;
			}
		}
	}
	
	/*! \brief unpack a grid
	 *
	 * \param mem preallocated memory from where to unpack the grid
	 * \param ps unpack-stat info
	 */
	template<int ... prp> inline void unpack(ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
		//if all of the aggregate properties are simple (don't have "pack()" member)
		if (has_pack_agg<T,prp...>::result::value == false)
		{
			//Call an unpacker
			unpack_simple_cond<sizeof...(prp) == 0,prp...>::unpack(*this,mem,ps);
		}
		//If at least one is not simple (has a "pack()" member)
		else
		{
			size_t dims[dim];
			
			//Unpack a size of a source grid
			for (size_t i = 0; i < dim; i++)
			{
				size_t u2 = 0;
				Unpacker<size_t, S>::unpack(mem,u2,ps);
				dims[i] = u2;
			}

			//Resize a destination grid
			this->resize(dims);
			//this->setMemory();

			auto key_it = this->getIterator();
	
			while (key_it.isNext())
			{
				auto k = key_it.get();
				//Call a pack request
				call_aggregateUnpack<decltype(this->get_o(k)),S,prp ... >::call_unpack(this->get_o(k),mem,ps);
				
				++key_it;
			}
		}
	}
	
	/*! \brief Pack finalize Finalize the pack of this object. In this case it does nothing
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void packFinalize(ExtPreAlloc<S> & mem, Pack_stat & sts)
	{}
	
	/*! \brief Pack the object into the memory given an iterator
	 *
	 * \tparam prp properties to pack
	 *
	 * \param mem preallocated memory where to pack the objects
	 * \param sub_it sub grid iterator ( or the elements in the grid to pack )
	 * \param sts pack statistic
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<S> & mem, grid_key_dx_iterator_sub<dims> & sub_it, Pack_stat & sts)
	{
#ifdef SE_CLASS1
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

		// Sending property object
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>,typename layout_base<prp_object>::type, layout_base,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(sub_it.getVolume());

		// copy all the object in the send buffer
		typedef encapc<dims,value_type,layout > encap_src;
		// destination object type
		typedef encapc<1,prp_object,typename dtype::layout_type > encap_dst;

		pack_with_iterator<sizeof...(prp) != T::max_prop || has_pack_gen<prp_object>::value,
						   dims,
						   decltype(*this),
						   encap_src,
						   encap_dst,
						   typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,
						   grid_key_dx_iterator_sub<dims>,
						   dtype,
						   prp...>::pack(*this,sub_it,dest);

		// Update statistic
		sts.incReq();
	}


	/*! \brief Insert an allocation request
	 *
	 * \tparam prp set of properties to pack
	 *

	 * \param sub sub-grid iterator
	 * \param vector of requests
	 *
	 */
	template<int ... prp> void packRequest(grid_key_dx_iterator_sub<dims> & sub, size_t & req)
	{
		typedef openfpm::vector<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type,ExtPreAlloc<S>,layout,layout_base,openfpm::grow_policy_identity> dtype;
		dtype dvect;

		// Calculate the required memory for packing
		size_t alloc_ele = dvect.template calculateMem<prp...>(sub.getVolume(),0);

		req += alloc_ele;
	}


	/*! \brief unpack the sub-grid object
	 *
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<unsigned int ... prp,typename S2, typename context_type> void unpack(ExtPreAlloc<S2> & mem, grid_key_dx_iterator_sub<dims> & sub_it, Unpack_stat & ps,context_type & context)
	{
		// object that store the information in mem
		typedef object<typename object_creator<typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin ,openfpm::grow_policy_identity> stype;

		size_t size = stype::template calculateMem(sub_it.getVolume(),0);

		// Create an object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

		// Create an object of the packed information over a pointer (No allocation is produced)
		stype src;
		src.setMemory(ptr);
		src.resize(sub_it.getVolume());

		// copy all the object in the send buffer
		typedef encapc<dims,grid_base_impl<dim,T,S,layout,layout_base>::value_type,layout > encap_dst;
		// destination object type
		typedef encapc<1,prp_object,typename memory_traits_lin<prp_object>::type > encap_src;

		unpack_with_iterator<dims,
							 decltype(*this),
							 encap_src,
							 encap_dst,
							 typename grid_base_impl<dim,T,S,layout,layout_base>::value_type::type,
							 grid_key_dx_iterator_sub<dims>,
							 stype,
							 prp...>::unpack(*this,sub_it,src);

		ps.addOffset(size);
	}

	/*! \brief unpack the sub-grid object applying an operation
	 *
	 * \tparam op operation
	 * \tparam prp properties to unpack
	 *
	 * \param mem preallocated memory from where to unpack the object
	 * \param sub sub-grid iterator
	 * \param obj object where to unpack
	 *
	 */
	template<template<typename,typename> class op, typename S2, unsigned int ... prp>
	void unpack_with_op(ExtPreAlloc<S2> & mem, grid_key_dx_iterator_sub<dim> & sub2, Unpack_stat & ps)
	{
		PtrMemory * ptr1;

		size_t sz[dim];

		for (size_t i = 0 ; i < dim ; i++)
			sz[i] = sub2.getStop().get(i) - sub2.getStart().get(i) + 1;

		size_t tot = 1;

		for (size_t i = 0 ; i < dim ; i++)
		{tot *= sz[i];}

		tot *= sizeof(T);

#ifdef SE_CLASS1

		if (ps.getOffset() + tot > mem.size())
			std::cerr << __FILE__ << ":" << __LINE__ << " Error: overflow in the receiving buffer for ghost_put" << std::endl;

#endif

		// add the received particles to the vector
		ptr1 = new PtrMemory(((char *)mem.getPointerBase()+ps.getOffset()),tot);

		// create vector representation to a piece of memory already allocated
		grid_base_impl<dim,T,PtrMemory,typename memory_traits_lin<T>::type,memory_traits_lin> gs;

		gs.setMemory(*ptr1);

		// resize with the number of elements
		gs.resize(sz);

		// Merge the information

		auto it_src = gs.getIterator();

		while (sub2.isNext())
		{
			object_s_di_op<op,decltype(gs.get_o(it_src.get())),decltype(this->get_o(sub2.get())),OBJ_ENCAP,prp...>(gs.get_o(it_src.get()),this->get_o(sub2.get()));

			++sub2;
			++it_src;
		}

		ps.addOffset(tot);
	}

	/*! \brief Calculate the memory size required to pack n elements
	 *
	 * Calculate the total size required to store n-elements in a vector
	 *
	 * \param n number of elements
	 * \param e unused
	 *
	 * \return the size of the allocation number e
	 *
	 */
	template<int ... prp> static inline size_t packMem(size_t n, size_t e)
	{
		if (sizeof...(prp) == 0)
			return n * sizeof(typename T::type);

		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		return n * sizeof(prp_object);
	}
	

