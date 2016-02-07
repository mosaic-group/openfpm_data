/// <summary>
/// SimpleRNG is a simple random number generator based on
/// George Marsaglia's MWC (multiply with carry) generator.
/// Although it is very simple, it passes Marsaglia's DIEHARD
/// series of random number generator tests.
///
/// Written by John D. Cook
/// http://www.johndcook.com
/// Converted to C++ by Pietro Incardona
/// </summary>

#include <string>
#include <cmath>

class SimpleRNG
{
private:
	uint m_w;
	uint m_z;

	// This is the heart of the generator.
	// It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
	// See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt
	uint GetUint()
	{
		m_z = 36969 * (m_z & 65535) + (m_z >> 16);
		m_w = 18000 * (m_w & 65535) + (m_w >> 16);
		return (m_z << 16) + m_w;
	}

public:

	SimpleRNG()
	{
		// These values are not magical, just the default values Marsaglia used.
		// Any pair of unsigned integers should be fine.
		m_w = 521288629;
		m_z = 362436069;
	}

	// The random generator seed can be set three ways:
	// 1) specifying two non-zero unsigned integers
	// 2) specifying one non-zero unsigned integer and taking a default value for the second
	// 3) setting the seed from the system time

	void SetSeed(uint u, uint v)
	{
		if (u != 0) m_w = u;
		if (v != 0) m_z = v;
	}

	void SetSeed(uint u)
	{
		m_w = u;
	}

	void SetSeedFromSystemTime()
	{
		long x = clock();
		SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
	}

	// Produce a uniform random sample from the open interval (0, 1).
	// The method will not return either end point.
	double GetUniform()
	{
		// 0 <= u < 2^32
		uint u = GetUint();
		// The magic number below is 1/(2^32 + 2).
		// The result is strictly between 0 and 1.
		return (u + 1.0) * 2.328306435454494e-10;
	}

	// Get normal (Gaussian) random sample with mean 0 and standard deviation 1
	double GetNormal()
	{
		// Use Box-Muller algorithm
		double u1 = GetUniform();
		double u2 = GetUniform();
		double r = std::sqrt( -2.0*std::log(u1) );
		double theta = 2.0*3.14159265359*u2;
		return r*std::sin(theta);
	}

	// Get normal (Gaussian) random sample with specified mean and standard deviation
	double GetNormal(double mean, double standardDeviation)
	{
		if (standardDeviation <= 0.0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Shape must be positive. Received." << std::to_string(standardDeviation);
			throw 100001;
		}
		return mean + standardDeviation*GetNormal();
	}

	// Get exponential random sample with mean 1
	double GetExponential()
	{
		return -std::log( GetUniform() );
	}

	// Get exponential random sample with specified mean
	double GetExponential(double mean)
	{
		if (mean <= 0.0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << "Mean must be positive. Received " << mean;
			throw 100002;
		}
		return mean*GetExponential();
	}

	double GetGamma(double shape, double scale)
	{
		// Implementation based on "A Simple Method for Generating Gamma Variables"
		// by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
		// Vol 26, No 3, September 2000, pages 363-372.

		double d, c, x, xsquared, v, u;

		if (shape >= 1.0)
		{
			d = shape - 1.0/3.0;
			c = 1.0/std::sqrt(9.0*d);
			for (;;)
			{
				do
				{
					x = GetNormal();
					v = 1.0 + c*x;
				}
				while (v <= 0.0);
				v = v*v*v;
				u = GetUniform();
				xsquared = x*x;
				if (u < 1.0 -.0331*xsquared*xsquared || std::log(u) < 0.5*xsquared + d*(1.0 - v + std::log(v)))
					return scale*d*v;
			}
		}
		else if (shape <= 0.0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Shape must be positive. Received" << shape << "\n";
			throw 100003;
		}
		else
		{
			double g = GetGamma(shape+1.0, 1.0);
			double w = GetUniform();
			return scale*g*std::pow(w, 1.0/shape);
		}
	}

	double GetChiSquare(double degreesOfFreedom)
	{
		// A chi squared distribution with n degrees of freedom
		// is a gamma distribution with shape n/2 and scale 2.
		return GetGamma(0.5 * degreesOfFreedom, 2.0);
	}

	double GetInverseGamma(double shape, double scale)
	{
		// If X is gamma(shape, scale) then
		// 1/Y is inverse gamma(shape, 1/scale)
		return 1.0 / GetGamma(shape, 1.0 / scale);
	}

	double GetWeibull(double shape, double scale)
	{
		if (shape <= 0.0 || scale <= 0.0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Shape and scale parameters must be positive. Recieved " << shape << " and " << scale;
			throw 100004;
		}
		return scale * std::pow(-std::log(GetUniform()), 1.0 / shape);
	}

	double GetCauchy(double median, double scale)
	{
		if (scale <= 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << "Scale must be positive. Received " << scale << "\n";
			throw 100005;
		}

		double p = GetUniform();

		// Apply inverse of the Cauchy distribution function to a uniform
		return median + scale*std::tan(3.14159265359*(p - 0.5));
	}

	double GetStudentT(double degreesOfFreedom)
	{
		if (degreesOfFreedom <= 0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Degrees of freedom must be positive. Received " << degreesOfFreedom;
			throw 100006;
		}

		// See Seminumerical Algorithms by Knuth
		double y1 = GetNormal();
		double y2 = GetChiSquare(degreesOfFreedom);
		return y1 / std::sqrt(y2 / degreesOfFreedom);
	}

	// The Laplace distribution is also known as the double exponential distribution.
	double GetLaplace(double mean, double scale)
	{
		double u = GetUniform();
		return (u < 0.5) ?
			mean + scale*std::log(2.0*u) :
			mean - scale*std::log(2*(1-u));
	}

	double GetLogNormal(double mu, double sigma)
	{
		return std::exp(GetNormal(mu, sigma));
	}

	double GetBeta(double a, double b)
	{
		if (a <= 0.0 || b <= 0.0)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << "Beta parameters must be positive. Received " << a << " and " << b << "\n";
			throw 100007;
		}

		// There are more efficient methods for generating beta samples.
		// However such methods are a little more efficient and much more complicated.
		// For an explanation of why the following method works, see
		// http://www.johndcook.com/distribution_chart.html#gamma_beta

		double u = GetGamma(a, 1.0);
		double v = GetGamma(b, 1.0);
		return u / (u + v);
	}
};
