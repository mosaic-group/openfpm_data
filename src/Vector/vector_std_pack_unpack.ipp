/*!
 * This file include the implemetation of packer and unpacker for std vector
 */


// Structures that do a nested packing, depending on the existence of 'pack' function inside of the object

	//There is no pack() inside
	template<bool cond, typename T1, typename Memory1, int ... prp>
	struct pack_cond
	{
		void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
		{
#ifdef DEBUG
			std::cout << "No pack() function inside! (packer)" << std::endl;
#endif
			//Call an array primitive packer
			Packer<openfpm::vector<T1>, Memory1, PACKER_ARRAY_PRIMITIVE>::pack(mem,obj,sts,obj.size());		
		}
	};

	//There is pack() inside
	template<typename T1, typename Memory1, int ... prp>
	struct pack_cond<true, T1, Memory1, prp...>
	{
		void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
		{
#ifdef DEBUG
			std::cout << "There is pack() function inside! (true, map_vector_std)" << std::endl;
#endif

			//Pack the size of a vector
			Packer<size_t, Memory1>::pack(mem,obj.size(),sts);

			//Call a packer in nested way
			for (int i = 0; i < obj.size(); i++) {
				obj.get(i).template pack<prp...>(mem,sts);
			}
		}
	};

	// Structures that do a nested unpacking, depending on the existence of 'pack' function inside of the object

	//There is no pack() inside
	template<bool cond, typename T1, typename Memory1, int ... prp>
	struct unpack_cond
	{
		void unpacking(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Unpack_stat & ps)
		{
#ifdef DEBUG
			std::cout << "No pack() function inside! (unpacker_std)" << std::endl;
#endif

			//Call the array of primitives unpacker
			Unpacker<openfpm::vector<T1>, Memory1, PACKER_ARRAY_PRIMITIVE>::unpack(mem,obj,ps);
		}

	};

	//There is pack() inside
	template<typename T1, typename Memory1, int ... prp>
	struct unpack_cond<true, T1, Memory1, prp...>
	{
		void unpacking(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Unpack_stat & ps)
		{
#ifdef DEBUG
			std::cout << "There is pack() function inside! (unpacker_std)" << std::endl;
#endif

			//Unpacking a size of a source vector
			size_t u1 = 0;
			Unpacker<size_t, Memory1>::unpack(mem,u1,ps);

			//Resize a destination vector
			obj.resize(u1);

			//Call an unpacker in nested way

			//std::cout<< demangle(typeid(obj.get(1)).name()) << std::endl;
			for (int i = 0; i < obj.size(); i++) {
				obj.get(i).template unpack<prp...>(mem,ps);
			}
		}
	};

	// Structures that do a pack request, depending on the existence of 'packRequest' function inside of the object

	//There is no packRequest() inside
	template<bool cond, typename T1, int ... prp>
	struct packRequest_cond
	{
		void packingRequest(openfpm::vector<T1> & obj, std::vector<size_t> & v)
		{
#ifdef DEBUG
			std::cout << "There is no packRequest() function inside! (packingRequest)" << std::endl;
#endif
				//openfpm::vector<T1>::calculateMem(obj.size(),0);
				//Pushback a size of number of elements of the internal vectors
				v.push_back(sizeof(obj.size()));

				size_t alloc_ele = obj.calculateMem<prp...>(obj.size(),0);

				v.push_back(alloc_ele);
		}

	};

	//There is packRequest() inside
	template<typename T1, int ... prp>
	struct packRequest_cond<true, T1, prp...>
	{
		void packingRequest(openfpm::vector<T1> & obj, std::vector<size_t> & v)
		{
#ifdef DEBUG
			std::cout << "There is packRequest() function inside! (packingRequest)" << std::endl;
#endif
			//Pushback a size of number of elements of the external vectors
			v.push_back(sizeof(obj.size()));

			//Call an packRequest in nested way
			for (int i = 0; i < obj.size(); i++)
			{
				obj.get(i).template packRequest<prp...>(v);
			}
		}
	};

	// Structures that calculate memory for an object, depending on the existence of 'calculateMem' function inside of the object

	//There is no calculateMem() inside
	template<bool cond, typename T1, int ... prp>
	struct calculateMem_cond
	{
		size_t calculateMemory(T1 & obj, size_t n, size_t e)
		{
#ifdef DEBUG
			std::cout << "There is no calculateMem() function inside! (calculateMemory)" << std::endl;
#endif
			return obj.size() * sizeof(T);
		}
	};

	//There is calculateMem() inside
	template<typename T1, int ... prp>
	struct calculateMem_cond<true, T1, prp...>
	{
		size_t calculateMemory(T1 & obj, size_t n, size_t e)
		{
#ifdef DEBUG
			std::cout << "There is calculateMem() function inside! (calculateMemory)" << std::endl;
#endif
			size_t res = 0;

			for (int i = 0; i < n; i++) {
				res += obj.get(i).calculateMem<prp...>(obj.get(i).size(),0);
			}

			return res;
		}
	};

	/*! \brief pack a vector selecting the properties to pack
	 *
	 * \param mem preallocated memory where to pack the vector
	 * \param obj object to pack
	 * \param sts pack-stat info
	 *
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

	static bool calculateMem()
	{
		return false;
	}


	template<int ... prp> void pack(ExtPreAlloc<HeapMemory> & mem, Pack_stat & sts)
	{
#ifdef DEBUG
		if (mem.ref() == 0)
			std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

#ifdef DEBUG
		std::cout << "Inside pack() function! (map_vector_std)" << std::endl;
#endif
		pack_cond<has_pack<T>::type::value, T, HeapMemory, prp...> p;
		p.packing(mem, *this, sts);

	}

	/*! \brief Insert an allocation request into the vector
	 *
	 * \param obj vector object to pack
	 * \param requests vector
	 *
	 */
	template<int ... prp> void packRequest(std::vector<size_t> & v)
	{
		packRequest_cond<has_packRequest<T>::value, T, prp...> pr;
		pr.packingRequest(*this, v);

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
	template<unsigned int ... prp> void unpack(ExtPreAlloc<HeapMemory> & mem, Unpack_stat & ps)
	{
#ifdef DEBUG
		std::cout << "Inside unpack function (std)" << std::endl;
#endif
		unpack_cond<has_pack<T>::type::value, T, HeapMemory, prp...> unp;
		unp.unpacking(mem, *this, ps);

	}
