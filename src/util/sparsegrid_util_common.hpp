#ifndef SPARSEGRID_UTIL_COMMON
#define SPARSEGRID_UTIL_COMMON

template<typename base>
struct sparse_grid_bck_wrapper_impl
{
    base & b;

    typedef typename base::value_type bt_type;

    sparse_grid_bck_wrapper_impl(base & b)
    :b(b)
    {}

    auto operator[](int i) -> decltype(b[0])
    {
        b[0];
    }

    template<typename T>
    sparse_grid_bck_wrapper_impl<base> & operator=(T c)
    {
        for (int i = 0; i < b.size() ; i++)
        {
            b[i] = c;
        }

        return *this;
    }

    operator bt_type() const
    {
        return b[0];
    }
};

// openfpm::detail::multi_array::sub_array_openfpm<T,N1,vmpl>

template<typename T,typename vmpl>
struct sparse_grid_bck_wrapper_impl<openfpm::detail::multi_array::sub_array_openfpm<T, 1, vmpl>>
{
    openfpm::detail::multi_array::sub_array_openfpm<T, 1, vmpl > b;

    
    sparse_grid_bck_wrapper_impl(openfpm::detail::multi_array::sub_array_openfpm<T, 1, vmpl > b)
    :b(b)
    {}

    auto operator[](int i) -> decltype(b[i][0])
    {
        return b[i][0];
    }

    auto operator[](int i) const -> decltype(b[i][0])
    {
        return b[i][0];
    }
};

template<typename chunk_def>
struct sparse_grid_bck_value
{
	chunk_def bck;

	sparse_grid_bck_value(chunk_def bck)
	:bck(bck)
	{}

	template<unsigned int p>
	auto get() -> sparse_grid_bck_wrapper_impl<typename std::remove_reference<decltype(bck.template get<p>())>::type>
	{
		return sparse_grid_bck_wrapper_impl<typename std::remove_reference<decltype(bck.template get<p>())>::type>(bck.template get<p>());
	}

	template<unsigned int p>
	auto get() const -> sparse_grid_bck_wrapper_impl<typename std::remove_reference<decltype(bck.template get<p>())>::type>
	{
		return sparse_grid_bck_wrapper_impl<typename std::remove_reference<decltype(bck.template get<p>())>::type>(bck.template get<p>());
	}
};

#endif
