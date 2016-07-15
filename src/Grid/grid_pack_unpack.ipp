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
	static inline void pack(const grid_cpu<dim,T,S,Mem> & obj, ExtPreAlloc<S> & mem, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif
		
		//Pack the size of a grid
		for (size_t i = 0; i < dim; i++)
		{
			Packer<size_t, S>::pack(mem,obj.getGrid().size(i),sts);
		}
		
#ifdef DEBUG
		std::cout << "Grid size is " << std::endl;	
		for (size_t i = 0; i < dim; i++)
		{
			std::cout << obj.getGrid().size(i) << std::endl;
		}
#endif
		
#ifdef DEBUG
		std::cout << "Inside pack_simple(not 0 prop) function! (grid_pack_unpack)" << std::endl;
#endif

		// Sending property object and vector
		typedef object<typename object_creator<typename grid_cpu<dim,T,S,Mem>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,ExtPreAlloc<S>,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());

		auto it = obj.getIterator();

		obj.pack_with_iterator<decltype(it),dtype,prp...>(it,dest);

		// Update statistic
		sts.incReq();
	}
};

//Without specified properties
template<int ... prp>
struct pack_simple_cond<true, prp ...>
{
	static inline void pack(const grid_cpu<dim,T,S,Mem> & obj, ExtPreAlloc<S> & mem, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif
		
		//Pack the size of a grid
		for (size_t i = 0; i < dim; i++)
		{
			Packer<size_t, S>::pack(mem,obj.getGrid().size(i),sts);
		}		
#ifdef DEBUG		
		std::cout << "Grid size is " << std::endl;	
		for (size_t i = 0; i < dim; i++)
		{
			std::cout << obj.getGrid().size(i) << std::endl;
		}
#endif		
		// Sending property object
		typedef openfpm::vector<T,ExtPreAlloc<S>,openfpm::grow_policy_identity> dtype;
		
#ifdef DEBUG
		std::cout << "Inside pack_simple(0 prop) function! (grid_pack_unpack)" << std::endl;
#endif
		
		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());
	
		auto obj_it = obj.getIterator();
	
		size_t id = 0;
		
		while (obj_it.isNext())
		{
			// Copy
			dest.get(id) = obj.get_o(obj_it.get());
			
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
	static inline void unpack(grid_cpu<dim,T,S,Mem> & obj , ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
		
#ifdef DEBUG
		std::cout << "Inside unpack_simple(not 0 prop) function! (grid_pack_unpack)" << std::endl;
#endif
	
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
		typedef object<typename object_creator<typename grid_cpu<dim,T,S,Mem>::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> stype;

		// Calculate the size to pack the object
		size_t size = obj.packMem<prp...>(obj.size(),0);
		
		// Create an object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

		// Create an object over a pointer (No allocation is produced)
		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());

		auto it = obj.getIterator();

		obj.unpack_with_iterator<decltype(it),stype,prp...>(it,src);

		ps.addOffset(size);
	}
};

//Without specified properties
template<int ... prp>
struct unpack_simple_cond<true, prp ...>
{
	static inline void unpack(grid_cpu<dim,T,S,Mem> & obj , ExtPreAlloc<S> & mem, Unpack_stat & ps)
	{
#ifdef DEBUG
		std::cout << "Inside unpack_simple(0 prop) function! (grid_pack_unpack)" << std::endl;
#endif
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
		typedef openfpm::vector<T,PtrMemory,openfpm::grow_policy_identity> stype;

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
			obj.get_o(obj_it.get()) = src.get(id);
		
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
		//std::cout << demangle(typeid(this).name()) << std::endl;
		//std::cout << this->size() << std::endl;
#ifdef DEBUG
		std::cout << "Inside grid_pack_unpack.ipp packRequest()" << std::endl;
#endif
		
		// If all of the aggregate properties do not have a "pack()" member	
		if (has_pack_agg<T,prp...>::result::value == false)
		{
#ifdef DEBUG
		std::cout << "All of the aggregate members are simple!(packRequest)" << std::endl;
#endif
			size_t alloc_ele = this->packMem<prp...>(this->size(),0);
			req += alloc_ele;
		}
		//If at least one property has "pack()"
		else
		{			
#ifdef DEBUG
			std::cout << "Not all of the aggregate members are simple!(packRequest)" << std::endl;
#endif
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
#ifdef DEBUG
		std::cout << "Inside grid_pack_unpack.ipp pack()" << std::endl;
#endif
		//If all of the aggregate properties are simple (don't have "pack()" member)
		if (has_pack_agg<T,prp...>::result::value == false)
		{
	#ifdef DEBUG
			std::cout << "All of the aggregate members are simple!(pack)" << std::endl;
	#endif
			//Call a packer
			pack_simple_cond<sizeof...(prp) == 0,prp...>::pack(*this,mem,sts);
		}
		//If at least one property has a "pack()" member
		else
		{
			
#ifdef DEBUG
			std::cout << "Not all of the aggregate members are simple!(pack)" << std::endl;
#endif
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
#ifdef DEBUG
			std::cout << "All of the aggregate members are simple!(unpack)" << std::endl;
#endif
			//Call an unpacker
			unpack_simple_cond<sizeof...(prp) == 0,prp...>::unpack(*this,mem,ps);
		}
		//If at least one is not simple (has a "pack()" member)
		else
		{
#ifdef DEBUG
			std::cout << "Not all of the aggregate members are simple!(unpack)" << std::endl;
#endif
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
	
	/*! \brief Pack an N-dimensional grid into a vector like structure B given an iterator of the grid
	 *
	 * \tparam it type of iterator of the grid-structure
	 * \tparam dtype type of the structure B
	 * \tparam dim Dimensionality of the grid
	 * \tparam properties to pack
	 *
	 * \param it Grid iterator
	 * \param obj object to pack
	 * \param dest where to pack
	 *
	 */
	template <typename it, typename dtype, int ... prp> void pack_with_iterator(it & sub_it, dtype & dest) const
	{
		// Sending property object
		typedef object<typename object_creator<typename grid_cpu<dim,T,S,Mem>::value_type::type,prp...>::type> prp_object;

		size_t id = 0;

		// Packing the information
		while (sub_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<dims,value_type,memory_conf > encap_src;
			// destination object type
			typedef encapc<1,prp_object,typename dtype::memory_conf > encap_dst;

			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(this->get_o(sub_it.get()),dest.get(id));

			++id;
			++sub_it;
		}
	}
	
	/*! \brief unpack the grid given an iterator
	 *
	 * \tparam it type of iterator
	 * \tparam prp of the grid object to unpack
	 *
	 */
	template <typename it, typename stype, unsigned int ... prp> void unpack_with_iterator(it & sub_it, stype & src)
	{
		size_t id = 0;

		// Sending property object
		typedef object<typename object_creator<typename grid_cpu<dim,T,S,Mem>::value_type::type,prp...>::type> prp_object;

		// unpacking the information
		while (sub_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<dims,grid_cpu<dim,T,S,Mem>::value_type,grid_cpu<dim,T,S,Mem>::memory_conf > encap_dst;
			// destination object type
			typedef encapc<1,prp_object,typename grid_cpu<dims,prp_object>::memory_conf > encap_src;

			// Copy only the selected properties
			object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),this->get_o(sub_it.get()));

			++id;
			++sub_it;
		}
	}
