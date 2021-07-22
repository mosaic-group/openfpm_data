/*!
 * This file contains the implemetation of packer and unpacker for vector
 * Created on: Jan 5, 2016
 *     Author: Yaroslav Zaluzhnyi and Pietro Incardona
 */

/*! This Function to indicate the vector class has a packer function
 *
 * \return true vector has a pack function
 *
 */
static bool pack()
{
	return true;
}

/*! This Function indicate that vector class has a packRequest function
 *
 * \return true vector has a packRequest function
 *
 */
static bool packRequest()
{
	return true;
}

/*! This Function indicate that vector class has a packMem function
 *
 * \return true vector has a packMem function
 *
 */
static bool packMem()
{
       return true;
}

//! Structures that do a nested packing, depending on the existence of 'pack()' function inside the object
//There is no pack() inside
template<bool cond, typename T1, typename Memory1, int ... prp>
struct pack_cond
{
	//! Trivial case it serialize the vector object as simply an array of objects
	void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
	{
		Packer<openfpm::vector<T1>, Memory1, PACKER_ARRAY_PRIMITIVE>::pack(mem,obj,sts,obj.size());	
	}

};

//! Structures that do a nested packing, depending on the existence of 'pack()' function inside the object
//There is pack() inside
template<typename T1, typename Memory1, int ... prp>
struct pack_cond<true, T1, Memory1, prp...>
{
	//! It traverse the type-tree structure to serialize the vector object into memory
	void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
	{
	   for (size_t i = 0; i < obj.size(); i++)
		   obj.get(i).template pack<prp...>(mem, sts);
	}
};

//! Structures that calculate how many bytes are required to serialize an object
// depending on the existence of 'packMem()' function inside of the object
// There is no packMem() inside
template<bool cond, typename T1>
struct packMem_cond
{
	//! calculate how many bytes are needed to serialize the object
	size_t packMemory(T1 & obj, size_t n, size_t e)
	{
		return grow_policy_double::grow(0,n) * sizeof(T);
	}

};

//! Structures that calculate memory for an object, depending on the existence of 'packMem()'
// function inside of the object
// There is packMem() inside
template<typename T1>
struct packMem_cond<true, T1>
{
	//! calculate how many bytes are needed to serialize the object
	size_t packMemory(T1 & obj, size_t n, size_t e)
	{
		size_t res = 0;
		size_t count = grow_policy_double::grow(0,n) - obj.size();
		for (size_t i = 0; i < n; i++) {
			res += obj.get(i).packMem(n,0);
		}
		return res+count*sizeof(T);
	}
};


//! These structures serialize a simple (no "pack()" inside) object
// With specified properties
template<bool sel, int ... prp>
struct pack_simple_cond
{
	//! serialize the vector
	template<typename Memory2> static inline void pack(const openfpm::vector<T,Memory,layout_base,grow_p,OPENFPM_NATIVE> & obj, ExtPreAlloc<Memory2> & mem, Pack_stat & sts)
	{
		//Pack the size of a vector
		Packer<size_t, Memory2>::pack(mem,obj.size(),sts);
		
		// Sending property object
		typedef openfpm::vector<T,Memory,layout_base,grow_p> vctr;
		typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;
	
		typedef openfpm::vector<prp_object,ExtPreAlloc<Memory2>, layout_base ,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());
	
		auto obj_it = obj.getIterator();
	
		while (obj_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<1,typename vctr::value_type,typename vctr::layout_type > encap_src;
			// destination object type
			typedef encapc<1,prp_object,typename dtype::layout_type > encap_dst;
	
			// Copy only the selected properties
			object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(obj.get(obj_it.get()),dest.get(obj_it.get()));
	
			++obj_it;
		}
	
		// Update statistic
		sts.incReq();	
	}
};

//! These structures serialize a simple (no "pack()" inside) object
// Without specified properties
template<int ... prp>
struct pack_simple_cond<true, prp ...>
{
	//! serialize the vector
	template<typename Memory2>
	static inline void pack(const openfpm::vector<T,Memory,layout_base,grow_p,OPENFPM_NATIVE> & obj , ExtPreAlloc<Memory2> & mem, Pack_stat & sts)
	{
		//Pack the size of a vector
		Packer<size_t, Memory2>::pack(mem,obj.size(),sts);
		
		// Sending property object
		typedef openfpm::vector<T,ExtPreAlloc<Memory2>,memory_traits_lin,openfpm::grow_policy_identity> dtype;

		// Create an object over the preallocated memory (No allocation is produced)
		dtype dest;
		dest.setMemory(mem);
		dest.resize(obj.size());
	
		auto obj_it = obj.getIterator();
	
		while (obj_it.isNext())
		{
			// Copy
			dest.get(obj_it.get()) = obj.get(obj_it.get());
	
			++obj_it;
		}
	
		// Update statistic
		sts.incReq();
	}
};

//! These structures do an de-serialize a simple object (no pack() inside)
// With specified properties
template<bool sel, int ... prp>
struct unpack_simple_cond
{
	//! De-serialize an object
	template<typename Memory2>
	static inline void unpack(openfpm::vector<T,Memory, layout_base,grow_p,OPENFPM_NATIVE> & obj , ExtPreAlloc<Memory2> & mem, Unpack_stat & ps)
	{
		//Unpack a size of a source vector
		size_t u2 = 0;
		Unpacker<size_t, Memory2>::unpack(mem,u2,ps);
		
		//Resize a destination vector
		obj.resize(u2);
		
		size_t id = 0;
		
		// Sending property object
		typedef openfpm::vector<T> vctr;
		typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;
		typedef openfpm::vector<prp_object,PtrMemory, memory_traits_lin,openfpm::grow_policy_identity> stype;


		// Calculate the size to pack the object
		size_t size = obj.packMem<prp...>(obj.size(),0);
		
		// Create a Pointer object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));
		
		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());
		auto obj_it = obj.getIterator();
		
		while (obj_it.isNext())
		{
			// copy all the object in the send buffer
			typedef encapc<1,typename vctr::value_type,typename vctr::layout_type > encap_dst;
			// destination object type
			typedef encapc<1,prp_object,typename stype::layout_type > encap_src;
		
			// Copy only the selected properties
			object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),obj.get(obj_it.get()));
		
			++id;
			++obj_it;
		}
		
		ps.addOffset(size);
	}
};

//! These structures de-serialize a simple object (no pack() inside)
//! unpack Without specified properties
template<int ... prp>
struct unpack_simple_cond<true, prp ...>
{
	/*! \brief unpack from the memory the data structure and put it into obj
	 *
	 * \param obj object to deserialize
	 * \param mem object containing the raw data to deserialize
	 * \param ps statistic
	 *
	 */
	template<typename Memory2>
	static inline void unpack(openfpm::vector<T,Memory,layout_base, grow_p,OPENFPM_NATIVE> & obj , ExtPreAlloc<Memory2> & mem, Unpack_stat & ps)
	{
		//Unpack a size of a source vector
		size_t u2 = 0;
		Unpacker<size_t, Memory2>::unpack(mem,u2,ps);
		
		//Resize a destination vector
		obj.resize(u2);
		
		size_t id = 0;
		
		// Sending property object
		typedef openfpm::vector<T,PtrMemory,memory_traits_lin,openfpm::grow_policy_identity> stype;

		// Calculate the size to pack the object
		size_t size = obj.packMem<prp...>(obj.size(),0);
		
		// Create a Pointer object over the preallocated memory (No allocation is produced)
		PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));
		
		stype src;
		src.setMemory(ptr);
		src.resize(obj.size());
		auto obj_it = obj.getIterator();
		
		while (obj_it.isNext())
		{
			// Copy
			obj.get(obj_it.get()) = src.get(id);
		
			++id;
			++obj_it;
		}
		
		ps.addOffset(size);
	}
};


/*! \brief It calculate the number of byte required to serialize the object
 *
 * \tparam prp list of properties
 *
 * \param req reference to the total counter required to pack the information
 *
 */
template<int ... prp> inline void packRequest(size_t & req) const
{
	//Pushback a sizeof number of elements of the internal vectors
	req += sizeof(this->size());
	
	// If all of the aggregate properties do not have a "pack()" member	
	if (has_pack_agg<T,prp...>::result::value == false)
	{
		size_t alloc_ele = this->packMem<prp...>(this->size(),0);
		req += alloc_ele;
	}
	//If at least one property has "pack()"
	else
	{
		for (size_t i = 0 ; i < this->size() ; i++)
		{
			//Call a pack request
			call_aggregatePackRequest<decltype(this->get(i)),Memory,prp ... >::call_packRequest(this->get(i),req);
		}
	}
}


/*! \brief pack a vector selecting the properties to pack
 *
 * \param mem preallocated memory where to pack the vector
 * \param sts pack-stat info
 *
 */
template<int ... prp> inline void pack(ExtPreAlloc<HeapMemory> & mem, Pack_stat & sts) const
{
	//If all of the aggregate properties are simple (don't have "pack()" member)
	if (has_pack_agg<T,prp...>::result::value == false)
	//if (has_aggregatePack<T,prp ... >::has_pack() == false)
	{
		//Call a packer
		pack_simple_cond<sizeof...(prp) == 0,prp...>::template pack(*this,mem,sts);
	}
	//If at least one property has a "pack()" member
	else
	{
		//Pack the size of a vector
		Packer<size_t, HeapMemory>::pack(mem,this->size(),sts);
		
		for (size_t i = 0 ; i < this->size() ; i++)
		{
			//Call a packer in nested way
			call_aggregatePack<decltype(this->get(i)),HeapMemory,prp ... >::call_pack(this->get(i),mem,sts);
		}
	}
}

/*! \brief unpack a vector
 *
 * \param mem preallocated memory from where to unpack the vector
 * \param ps unpack-stat info
 */
template<int ... prp> inline void unpack(ExtPreAlloc<HeapMemory> & mem, Unpack_stat & ps)
{
	//if all of the aggregate properties are simple (don't have "pack()" member)
	if (has_pack_agg<T,prp...>::result::value == false)
	//if (has_aggregatePack<T,prp ... >::has_pack() == false)
	{
		//Call an unpacker
		unpack_simple_cond<sizeof...(prp) == 0,prp...>::unpack(*this,mem,ps);
	}
	//If at least one is not simple (has a "pack()" member)
	else
	{
		//Unpack a size of a source vector
		size_t u2 = 0;
		Unpacker<size_t, HeapMemory>::unpack(mem,u2,ps);
		
		//Resize a destination vector
		this->resize(u2);
		
		for (size_t i = 0 ; i < this->size() ; i++)
		{
			//Call an unpacker in nested way
			call_aggregateUnpack<decltype(this->get(i)),HeapMemory,prp ... >::call_unpack(this->get(i),mem,ps);
		}
	}
}

