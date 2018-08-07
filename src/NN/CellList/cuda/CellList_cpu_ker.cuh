/*
 * CellList_gpu_ker.cuh
 *
 *  Created on: Jul 30, 2018
 *      Author: i-bird
 */

#ifndef CELLLIST_CPU_KER_CUH_
#define CELLLIST_CPU_KER_CUH_

template<unsigned int dim, typename T, typename Mem_type, typename transform>
class CellList_cpu_ker: Mem_type
{

public:

	CellList_cpu_ker(const Mem_type & mt)
	:Mem_type(mt)
	{}

};


#endif /* CELLLIST_GPU_KER_CUH_ */
