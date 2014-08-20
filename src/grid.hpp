#include <iostream>
#include <boost/shared_array.hpp>
#include <vector>
#include <initializer_list>
#include <array>

#define HARDWARE 1

template<typename T, unsigned int s1, unsigned int s2, unsigned int s3>
class tensor
{
  T g_mem[s1][s2][s3];
 
public:
  tensor() {};
  ~tensor() {};
  
  size_t size()
  {
    return s1*s2*s3*g_mem[0][0][0].size();
  };
};


template<unsigned int dim>
class grid_key_dx
{
public:
  
  grid_key_dx()
  {
  }
  
  template<typename a, typename ...T>grid_key_dx(a v,T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  template<typename a, typename ...T>void set(a v, T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  mem_id k[dim];  
  
private:
  
  template<typename a, typename ...T>void invert_assign(a v,T...t)
  {
    k[sizeof...(T)] = v;
    invert_assign(t...);
  }

  template<typename a, typename ...T>void invert_assign(a v)
  {
    k[0] = v;
  }

};

template<unsigned int dim, unsigned int p>
class grid_key_d
{
public:
  
  
  template<typename a, typename ...T>grid_key_d(a v,T...t)
  {
    k[dim-1] = v;
    invert_assign(t...);
  }
  
  template<typename a, typename ...T>void invert_assign(a v,T...t)
  {
    k[sizeof...(T)] = v;
    invert_assign(t...);
  }

  template<typename a, typename ...T>void invert_assign(a v)
  {
    k[0] = v;
  }
  
  mem_id k[dim];
};

template<unsigned int p>
class grid_key_1
{
public:
  
  grid_key_1(int i_)
  {
    k[0] = i_;
  }
  
  mem_id k[1];
};

template<unsigned int p>
class grid_key_2
{
public:
  
  grid_key_2(int i_, int j_)
  {
    k[0] = j_;
    k[1] = i_;
  }
  
  mem_id k[2];
};

template<unsigned int p>
class grid_key_c3
{
public:
  
  grid_key_c3(mem_id i_, mem_id j_, mem_id k_)
  :k{{k_,j_,i_}}
  {
  }
  
  const std::array<mem_id,3> k;
};

template<unsigned int p>
class grid_key_3
{
public:
  
  grid_key_3(mem_id i_, mem_id j_, mem_id k_)
  {
    k[0] = k_;
    k[1] = j_;
    k[2] = i_;
  }
  
  mem_id k[3];
};

template<unsigned int p>
class grid_key_4
{
public:
  
  grid_key_4(mem_id u_, mem_id i_, mem_id j_, mem_id k_)
  {
    k[0] = k_;
    k[1] = j_;
    k[2] = i_;
    k[3] = u_;
  }
  
  mem_id k[4];
};

template<unsigned int p>
class grid_key
{
public:
  
  template<unsigned int s>grid_key(grid_key<s> & key)
  {
    id = key.id;
  }
  
  grid_key(size_t sz)
  :id(new mem_id[sz])
  {
  }
  
  void set_d(size_t id_, mem_id val)
  {
    id[id_] = val;
  }
  
  void set(mem_id val1)
  {
    id[0] = val1;
  }
  
  void set(mem_id val1, mem_id val2)
  {
    id[0] = val2;
    id[1] = val1;
  }
  
  void set(mem_id val1, mem_id val2, mem_id val3)
  {
    id[0] = val3;
    id[1] = val2;
    id[2] = val1;
  }
  
  void set(mem_id val1, mem_id val2, mem_id val3, mem_id val4)
  {
    id[0] = val4;
    id[1] = val3;
    id[2] = val2;
    id[3] = val1;
  }
  
  mem_id * getId()
  {
    return id.get();
  }
  
  boost::shared_array<mem_id> id;
};

//#pragma openfpm create(layout)
template<unsigned int N, typename T>
class grid
{
  size_t size_tot;
  size_t sz[N];
  size_t sz_s[N];
  
public:
  
  // Static element to calculate total size

  size_t totalSize(std::vector<size_t> & sz)
  {
    size_t tSz = 1;
    
    for (size_t i = 0 ;  i < sz.size() ; i++)
    {
      tSz *= sz[i];
    }
    
    return tSz;
  }
  
  grid(std::vector<size_t> & sz)
  : size_tot(totalSize(sz))
  {
    sz_s[0] = sz[0];
    this->sz[0] = sz[0];
    for (size_t i = 1 ;  i < sz.size() ; i++)
    {
      sz_s[i] = sz[i]*sz_s[i-1];
      this->sz[i] = sz[i];
    }
  }
  
  #pragma openfpm layout(get)
  template<unsigned int dim> mem_id LinId(grid_key_dx<dim> & gk)
  {
    mem_id lid = gk.k[0];
    for (mem_id i = 1 ; i < dim ; i++)
    {
      lid += gk.k[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  template<unsigned int p> mem_id LinId(grid_key_3<p> & gk)
  {
    mem_id lid = gk.k[0];
    lid += gk.k[1]*sz_s[0];
    lid += gk.k[2]*sz_s[1];
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  template<unsigned int dim, unsigned int p> mem_id LinId(grid_key_d<dim,p> & gk)
  {
    mem_id lid = gk.k[0];
    for (mem_id i = 1 ; i < dim ; i++)
    {
      lid += gk.k[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  mem_id LinId(mem_id * id)
  {
    mem_id lid = 0;
    lid += id[0];
    for (mem_id i = 1 ; i < N ; i++)
    {
      lid += id[i] * sz_s[i-1];
    }
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  mem_id LinId(mem_id i)
  {
    return i;
  }
  
  #pragma openfpm layout(get)
  mem_id LinId(mem_id j, mem_id k)
  {
    mem_id lid = 0;
    lid += k;
    lid += j*sz_s[0];
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  mem_id LinId(mem_id i, mem_id j, mem_id k)
  {
    mem_id lid = 0;
    lid += k;
    lid += j*sz_s[0];
    lid += i*sz_s[1];
    
    return lid;
  }
  
  #pragma openfpm layout(get)
  mem_id LinId(mem_id u, mem_id i, mem_id j, mem_id k)
  {
    mem_id lid = 0;
    lid += k;
    lid += j*sz_s[0];
    lid += i*sz_s[1];
    lid += u*sz_s[2];
    
    return lid;
  }
  
  ~grid() {};
  
  #pragma openfpm layout(size)
  size_t size()
  {
    return size_tot;
  };

  /**
   *
   * Get the size of the grid on direction i
   *
   * @param i direction
   *
   */

  size_t size(unsigned int i)
  {
	  return sz_s[i];
  }
};

/**
 *
 * Grid key class iterator, iterate through the grid_key
 *
 */

template<unsigned int dim>
class grid_key_dx_iterator
{
	grid<dim,void> & grid_base;
	template<typename T> grid_key_dx_iterator(grid<dim,T> & g)
	: grid_base(g)
	{}

	grid_key_dx<dim> gk;

	grid_key_dx<dim> next()
	{
		size_t id = gk.get(0);
		gk.set(0,id+1);
		for (int i = 0 ; i < dim ; i++)
		{
			size_t id = gk.get(i);
			if (id == grid_base.size(i))
			{
				gk.set(i,0);
				id = gk.get(i+1);
				gk.set(i+1,id+1);
			}
			else
			{
				break;
			}
		}
		return gk++;
	}

	bool hasNext()
	{
		for (int i = dim-1 ; i >= 0 ; i--)
		{
			if (gk.get(i) < grid_base.size(i))
			{
				return false;
			}
		}

		return true;
	}
};

template<unsigned int s1, unsigned int s2, unsigned int s3>
class tensor<int,s1,s2,s3>
{
  
public:
  size_t size()
  {
    return s1*s2*s3;
  }
};






