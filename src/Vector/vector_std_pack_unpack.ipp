/*!
 * This file contains the implemetation of packer and unpacker for std vector
 * Created on: Jan 5, 2016
 *     Author: Yaroslav Zaluzhnyi
 */


	//! Structures that do a nested packing, depending on the existence of 'pack' function inside of the object
	//There is no pack() inside
	template<bool cond, typename T1, typename Memory1, int ... prp>
	struct pack_cond
	{
		//! Serialize the object
		void packing(ExtPreAlloc<Memory1> & mem, const openfpm::vector<T1> & obj, Pack_stat & sts)
		{
			//Call an array primitive packer
			Packer<openfpm::vector<T1>, Memory1, PACKER_ARRAY_PRIMITIVE>::pack(mem,obj,sts,obj.size());		
		}
	};

	//! Structures that do a nested packing, depending on the existence of 'pack' function inside of the object
	//There is pack() inside
	template<typename T1, typename Memory1, int ... prp>
	struct pack_cond<true, T1, Memory1, prp...>
	{
		//! Serialize the object
		void packing(ExtPreAlloc<Memory1> & mem, const openfpm::vector<T1> & obj, Pack_stat & sts)
		{
			//Pack the size of a vector
			Packer<size_t, Memory1>::pack(mem,obj.size(),sts);

			//Call a packer in nested way
			for (size_t i = 0; i < obj.size(); i++) {
				obj.get(i).template pack<prp...>(mem,sts);
			}
		}
	};

	//! Structures that do a nested unpacking, depending on the existence of 'pack' function inside the object
	// There is no pack() inside
	template<bool cond, typename T1, typename Memory1, int ... prp>
	struct unpack_cond
	{
		//! De-serialize the object
		void unpacking(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Unpack_stat & ps)
		{
			//Call the array of primitives unpacker
			Unpacker<openfpm::vector<T1>, Memory1, PACKER_ARRAY_PRIMITIVE>::unpack(mem,obj,ps);
		}

	};

	//! Structures that do a nested unpacking, depending on the existence of 'pack' function inside the object
	// There is pack() inside
	template<typename T1, typename Memory1, int ... prp>
	struct unpack_cond<true, T1, Memory1, prp...>
	{
		//! De-serialize the object
		void unpacking(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Unpack_stat & ps)
		{
			//Unpacking a size of a source vector
			size_t u1 = 0;
			Unpacker<size_t, Memory1>::unpack(mem,u1,ps);

			//Resize a destination vector
			obj.resize(u1);

			//Call an unpacker in nested way

			//std::cout<< demangle(typeid(obj.get(1)).name()) << std::endl;
			for (size_t i = 0; i < obj.size(); i++) {
				obj.get(i).template unpack<prp...>(mem,ps);
			}
		}
	};

	//! Structures that do a pack request, depending on the existence of 'packRequest' function inside of the object
	//There is no packRequest() inside
	template<bool cond, typename T1, int ... prp>
	struct packRequest_cond
	{
		//! Calculate the size required to serialize the object
		void packingRequest(const openfpm::vector<T1> & obj, size_t & req)
		{
				//Pushback a size of number of elements of the internal vectors
				req += sizeof(obj.size());

				size_t alloc_ele = obj.template packMem<prp...>(obj.size(),0);

				req += alloc_ele;
		}

	};

	//! Structures that do a pack request, depending on the existence of 'packRequest' function inside of the object
	//There is packRequest() inside
	template<typename T1, int ... prp>
	struct packRequest_cond<true, T1, prp...>
	{
		//! Calculate the size required to serialize the object
		void packingRequest(const openfpm::vector<T1> & obj, size_t & req)
		{
			//Pushback a size of number of elements of the external vectors
			req += sizeof(obj.size());

			//Call an packRequest in nested way
			for (size_t i = 0; i < obj.size(); i++)
			{
				obj.get(i).template packRequest<prp...>(req);
			}
		}
	};

	//! Structures that calculate memory for an object, depending on the existence of 'packMem' function inside the object
	//There is no packMem() inside
	template<bool cond, typename T1, int ... prp>
	struct packMem_cond
	{
		//! Calculate the bytes required to serialize the vector
		size_t packMemory(const T1 & obj, size_t n, size_t e)
		{
			return obj.size() * sizeof(T);
		}
	};

	//! Structures that calculate memory for an object, depending on the existence of 'packMem' function inside the object
	// There is packMem() inside
	template<typename T1, int ... prp>
	struct packMem_cond<true, T1, prp...>
	{
		//! Calculate the bytes required to serialize the vector
		size_t packMemory(const T1 & obj, size_t n, size_t e)
		{
			size_t res = 0;

			for (size_t i = 0; i < n; i++) {
				res += obj.get(i).template packMem<prp...>(obj.get(i).size(),0);
			}

			return res;
		}
	};
	

	//! Indicate that the vector has a packer function
	static bool pack()
	{
		return true;
	}

	//! Indicate that the vector has the function packRequest
	static bool packRequest()
	{
		return true;
	}

	//! Indicate that the vector has the function pack Memory
	static bool packMem()
	{
	       return true;
	}

	
	/*! \brief pack a vector
	 *
	 * \param mem preallocated memory where to pack the vector
	 * \param sts pack-stat info
	 *
	 */
	template<int ... prp> void pack(ExtPreAlloc<HeapMemory> & mem, Pack_stat & sts) const
	{
		//Call a nested packer
		pack_cond<has_pack<T>::type::value, T, HeapMemory, prp...> p;
		p.packing(mem, *this, sts);

	}

	/*! \brief Insert an allocation request into the vector
	 *
	 * \tparam prp properties to pack
	 *
	 * \param req counter with the byte required to pack the vector
	 *
	 */
	template<int ... prp> void packRequest(size_t & req) const
	{
		//Call a nested pack request
		packRequest_cond<has_packRequest<T>::value, T, prp...> pr;
		pr.packingRequest(*this, req);

	}

	/*! \brief unpack a vector
	 *
	 * \warning the properties should match the packed properties,
	 *
	 * \param mem preallocated memory from where to unpack the vector
	 * \param ps unpack info
	 *
	 */
	template<unsigned int ... prp> void unpack(ExtPreAlloc<HeapMemory> & mem, Unpack_stat & ps)
	{
		unpack_cond<has_pack<T>::type::value, T, HeapMemory, prp...> unp;
		unp.unpacking(mem, *this, ps);

	}
	
	/*! \brief Save this object into file
	 * 
	 * \param file filename
	 * 
	 * \return true if succed
	 * 
	 */
	bool save(const std::string & file) const
	{
		size_t req = 0;
		
		Packer<openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR>,HeapMemory>::packRequest(*this,req);

		// allocate the memory
		HeapMemory pmem;
		//pmem.allocate(req);
		ExtPreAlloc<HeapMemory> mem(req,pmem);

		//Packing

		Pack_stat sts;
		Packer<openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR>,HeapMemory>::pack(mem,*this,sts);

		// Save into a binary file
	    std::ofstream dump (file, std::ios::out | std::ios::binary);
	    if (dump.is_open() == false)
	    	return false;
	    dump.write ((const char *)pmem.getPointer(), pmem.size());
	    
	    return true;
	}
	
	/*! \brief Load this object from file
	 * 
	 * \param file filename
	 * 
	 * \return true if succed
	 * 
	 */
	bool load(const std::string & file)
	{
	    std::ifstream fs (file, std::ios::in | std::ios::binary | std::ios::ate );
	    if (fs.is_open() == false)
	    	return false;
	    
	    // take the size of the file
	    size_t sz = fs.tellg();
	    
	    fs.close();
	    
	    // reopen the file without ios::ate to read
	    std::ifstream input (file, std::ios::in | std::ios::binary );
	    if (input.is_open() == false)
	    	return false;
	    
	    // Create the HeapMemory and the ExtPreAlloc memory
	    size_t req = 0;
	    req += sz;
	    HeapMemory pmem;
		ExtPreAlloc<HeapMemory> mem(req,pmem);
		
		// read
	    input.read((char *)pmem.getPointer(), sz);
	    
	    //close the file
	    input.close();
		
		//Unpacking
		Unpack_stat ps;

	 	Unpacker<openfpm::vector<T,HeapMemory,layout,memory_traits_lin,grow_policy_double,STD_VECTOR>,HeapMemory>::unpack(mem,*this,ps);
	 	
	 	return true;
	}
