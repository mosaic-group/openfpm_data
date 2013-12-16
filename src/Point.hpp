

template<typename T> class Point
{
  T x;
  T y;
  T z;
  
  inline void setx(T x_)	{x = x_;};
  inline void sety(T y_)	{y = y_;};
  inline void setz(T z_)	{z = z_;};
};