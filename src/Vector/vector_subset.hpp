#ifndef VECTOR_SUBSET_HPP
#define VECTOR_SUBSET_HPP

namespace openfpm
{

template<unsigned int dim,
         typename prop,
         template<typename> class layout_base = memory_traits_inte>
    class vector_subset_ker
    {
        mutable openfpm::vector_gpu_ker<typename apply_transform<layout_base,prop>::type,layout_base> v_all;

        mutable openfpm::vector_gpu_ker<aggregate<int>,layout_base> indexes;

        public:

        vector_subset_ker(openfpm::vector_gpu_ker<typename apply_transform<layout_base,prop>::type,layout_base> & v_all,
                          openfpm::vector_gpu_ker<aggregate<int>,layout_base> & indexes)
        :v_all(v_all),indexes(indexes)
        {}

        // get the

		template <unsigned int p>
		__device__ __host__ inline auto get(size_t id) const -> decltype(v_all.template get<p>(0))
		{
            return v_all.template get<p>(indexes.template get<0>(id));
        }
    }

public:

    template<typename T, typename Memory, template<typename> class layout_base, typename grow_p, unsigned int impl>
    class vector_sub
    {
        vector<T,Memory,layout_base,grow_p,impl> & v_all;

        vector<aggregate<int>,Memory,layout_base,grow_p,impl> & indexes;

        public:

            vector(vector<T,Memory,layout_base,grow_p,impl> & v_all, 
                   vector<aggregate<int>,Memory,layout_base,grow_p,impl> & indexes)
            :v_all(v_all),indexes(indexes)
            {}


            // To kernel
            template<unsigned int ... prp> vector_dist_ker<dim,St,prop,layout_base> toKernel()
            {
                vector_subset_ker<dim,St,prop,layout_base> v(v_all.toKernel(), indexes.toKernel());

                return v;
            }

		    template <unsigned int p>
		    inline auto get(size_t id) const -> decltype(v_all.template get<p>(0))
		    {
                return v_all.template get<p>(indexes.template get<0>(id));
            }

    }

};

#endif