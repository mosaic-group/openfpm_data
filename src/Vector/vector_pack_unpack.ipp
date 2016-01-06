/*!
 * This file include the implemetation of packer and unpacker for vector
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

static bool packMem()
{
	return false;
}



// Structures that do a nested packing, depending on the existence of 'pack' function inside of the object

//There is no pack() inside
template<bool cond, typename T1, typename Memory1, int ... prp>
struct pack_cond
{
	void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
	{
		//std::cout << "No pack() function inside!" << std::endl;
	}

};


//There is pack() inside
template<typename T1, typename Memory1, int ... prp>
struct pack_cond<true, T1, Memory1, prp...>
{
	void packing(ExtPreAlloc<Memory1> & mem, openfpm::vector<T1> & obj, Pack_stat & sts)
	{
	   for (int i = 0; i < obj.size(); i++)
		   obj.get(i).template pack<prp...>(mem, sts);
	}
};

// Structures that calculate memory for an object, depending on the existence of 'packMem' function inside of the object

//There is no packMem()
template<bool cond, typename T1>
struct packMem_cond
{
	size_t packMemory(T1 & obj, size_t n, size_t e)
	{
#ifdef DEBUG
		std::cout << grow_policy_double::grow(0,n) << "*" << sizeof(T) << " " << demangle(typeid(T).name()) << std::endl;
#endif
		return grow_policy_double::grow(0,n) * sizeof(T);
	}

};

//There is packMem()
template<typename T1>
struct packMem_cond<true, T1>
{
	size_t packMemory(T1 & obj, size_t n, size_t e)
	{
		size_t res = 0;
		size_t count = grow_policy_double::grow(0,n) - obj.size();
		for (int i = 0; i < n; i++) {
			res += obj.get(i).packMem(n,0);
		}
#ifdef DEBUG
		std::cout << "Inside packMem - true (map_vector)" << std::endl;
#endif
		return res+count*sizeof(T);
	}
};

/*! \brief pack a vector selecting the properties to pack
 *
 * \param mem preallocated memory where to pack the vector
 * \param obj object to pack
 * \param sts pack-stat info
 *
 */
template<int ... prp> void pack(ExtPreAlloc<Memory> & mem, Pack_stat & sts)
{
#ifdef DEBUG
	if (mem.ref() == 0)
		std::cerr << "Error : " << __FILE__ << ":" << __LINE__ << " the reference counter of mem should never be zero when packing \n";
#endif

	//Pack the size of a vector
	Packer<size_t, Memory>::pack(mem,this->size(),sts);
	// Sending property object
	typedef openfpm::vector<T> vctr;
	typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;

	typedef openfpm::vector<prp_object,ExtPreAlloc<Memory>,openfpm::grow_policy_identity> dtype;
#ifdef DEBUG
	std::cout << "Inside pack() function! (map_vector)" << std::endl;
#endif
	// Create an object over the preallocated memory (No allocation is produced)
	dtype dest;
	dest.setMemory(mem);
	dest.resize(this->size());

	auto obj_it = this->getIterator();

	while (obj_it.isNext())
	{
		// copy all the object in the send buffer
		typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_src;
		// destination object type
		typedef encapc<1,prp_object,typename dtype::memory_conf > encap_dst;

		// Copy only the selected properties
		object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(this->get(obj_it.get()),dest.get(obj_it.get()));

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
template<int ... prp> void packRequest(std::vector<size_t> & v)
{

	//Pushback a size of number of elements of the internal vectors
	v.push_back(sizeof(this->size()));
#ifdef DEBUG
	std::cout << "Inside map_vector.hpp packRequest()" << std::endl;
#endif

	size_t alloc_ele = this->packMem<prp...>(this->size(),0);

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
template<unsigned int ... prp> void unpack(ExtPreAlloc<Memory> & mem, Unpack_stat & ps)
{

	//Unpack a size of a source vector
	size_t u2 = 0;
	Unpacker<size_t, Memory>::unpack(mem,u2,ps);

	//Resize a destination vector
	this->resize(u2);

	size_t id = 0;

	// Sending property object
	typedef openfpm::vector<T> vctr;
	typedef object<typename object_creator<typename vctr::value_type::type,prp...>::type> prp_object;
	typedef openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> stype;
	stype svect;

#ifdef DEBUG
	std::cout << "Inside unpack() function! (map_vector)" << std::endl;
#endif
	// Calculate the size to pack the object
	size_t size = this->packMem<prp...>(this->size(),0);

	// Create a Pointer object over the preallocated memory (No allocation is produced)
	PtrMemory & ptr = *(new PtrMemory(mem.getPointerOffset(ps.getOffset()),size));

	stype src;
	src.setMemory(ptr);
	src.resize(this->size());
	auto obj_it = this->getIterator();

	while (obj_it.isNext())
	{
		// copy all the object in the send buffer
		typedef encapc<1,typename vctr::value_type,typename vctr::memory_conf > encap_dst;
		// destination object type
		typedef encapc<1,prp_object,typename stype::memory_conf > encap_src;

		// Copy only the selected properties
		object_s_di<encap_src,encap_dst,OBJ_ENCAP,prp...>(src.get(id),this->get(obj_it.get()));

		++id;
		++obj_it;
	}

	ps.addOffset(size);
}

