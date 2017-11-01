// standard_brownian.h - Standard  Brownian motion.
#pragma once
#include "option.h"

namespace option {
#pragma warning(disable: 4355)

template<class R = double, class T = double, class S = double, class Z = double>
struct standard_brownian : public distribution_base<R,T,S,Z> { // B_t
	standard_brownian(const T& t)
		: distribution_base(t)
	{ }
	~standard_brownian()
	{ }
	standard_brownian& operator=(const standard_brownian& t)
	{
		if (this != &t)
			t_ = t.t_;

		return *this;
	}

	// probability density function
	static R static_pdf(const Z& z)
	{
		return t_ ? normal_pdf(z/sqrt(t_))/sqrt(t_) : 0.; // really delta function at 0
	}
	R pdf(const Z& z) const
	{
		return static_pdf(z);
	}

	// cumulative distribution function
	static R static_cdf(const T& t, const Z& z)
	{
		return t != 0 ? normal_cdf(z/sqrt(t)) : 1*(z >= 0.);
	}
	R cdf(const Z& z) const
	{
		return static_cdf(t_, z);
	}

	// kappa(t, s) = log E[exp(s X_t)]
	static R static_cumulant(const T& t, const S& s)
	{
		return t*s*s/2;
	}
	R cumulant(const S& s) const
	{
		return static_cumulant(t_, s);
	}

	// cumulant distribution function
	// log E[exp(s B_t) 1(B_t <= z)] = s^2t/2 + log P(X_t + st <= z)
	static R static_cudf(const T& t, const S& s, const Z& z)
	{
		return cumulant(s) + log(cdf(z - s*t));
	}
	R cudf(const S& s, const Z& z) const
	{
		return static_cudf(t_, s, z);
	}

	// exp(-kappa(s) + kappa(s, z))
	static R static_N(const T& t, const S& s, const Z& z)
	{
		return static_cdf(t, z - s*t);
	}
	R N(const S& s, const Z& z) const
	{
		return static_N(t_, s, z);
	}

};

	template<class T>
	T implied_guess(const T& f, const T& p, const T& k, const T& t)
	{
		if (k > 0)
			return implied_guess(f, p - f + k, -k, t);

		T z = sqrt(t)*normal_inv(p/(-k - f));

		return (z + sqrt(z*z - 2*t*log(-k/f)))/t;
	}
	template<class R>
		template<class T>
		static R black<R>::implied(const R& f, const R& p, const R& k, const T& t)
	{
		R s0 = 0.2, eps = 1e-11;
		int max_iteration_count = 100;
		R c = 1, s1, p1;
		
		eps *= p;
		ensure (eps != 0);

		if (k < 0) {
			return implied(f, f - k + p, -k, t);
		}

		// ensure price in 0 - infty vol range
		ensure (p >= __max(c*(f - k),0.));
		ensure ((c == 1 && p < f) || (c == -1 && p < k));

		R p0 = black<R>::value(f, s0, c*k, t) - p;

		// lucky guess
		if (fabs(p0) < eps)
			return s0;

		// bracket the root
		R m = 1.4;
		if (p0 > 0) {
			s1 = s0/m;
			p1 = black<R>::value(f, s1, c*k, t) - p;
			while (p1 > 0) {
				s0 = s1;
				p0 = p1;
				s1 = s0/m;
				p1 = black<R>::value(f, s1, c*k, t) - p;
			}
		}
		else {
			s1 = s0*m;
			p1 = black<R>::value(f, s1, c*k, t) - p;
			while (p1 < 0) {
				s0 = s1;
				p0 = p1;
				s1 = s0*m;
				p1 = black<R>::value(f, s1, c*k, t) - p;
			}
		}

		if (fabs(p1) < eps)
			return s1;

		ensure (p0*p1 < 0);

		// polish
		R ds = 0;
		R s2 = s0 - p0*(s1 - s0)/(p1 - p0);
		R p2 = black<R>::value(f, s2, c*k, t) - p;
		ds = black<R>::vega(f, s2, c*k, t);

		// if sigma is too small use bisection
		if (ds < 1e-4) {
			while (fabs(p2) > eps) {
				if (p0*p2 < 0) {
					s1 = s2;
				}
				else {
					ensure (p1*p2 < 0);
					s0 = s2;
				}
				s2 = (s1 + s0)/2;
				p2 = black<R>::value(f, s2, c*k, t) - p;
			}

			return s2;
		}

		// Newton-Raphson
		s0 = s2;
		p0 = p2;
		for (int i = 0; fabs(p0) > eps; ++i) {
			ensure (i < max_iteration_count);
			ensure (ds != 0);
			
			s1 = s0 - p0/ds;
			if (s1 < 0)
				s1 = s0/2;
			ds = 0;
			s0 = s1;
			p0 = black<R>::value(f, s0, c*k, t) - p;
			ds = black<R>::vega(f, s2, c*k, t);

		}

		return s0;
	}
	
} // namespace option


