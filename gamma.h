// gamma.h - Gamma process
// X_1 ~ log Gamma(a, b)
// E[Gamma(a,b)] = ab
// Var Gamma(a,b) = ab^2
// F = f exp(-kappa(t,s) + s X_t), X_t ~ log Gamma(a,b)
#pragma once
#include "option.h"
#include "../fmsdual/gamma.h"

// Gamma(alpha, beta) ~ exp(Normal(0,1))
// Mean and variance match for these values

#define GAMMA_ALPHA 0.5819767068693264243850020051090115585468693010753961362667 // 1/(e - 1)
#define GAMMA_BETA  0.9595173756674718597461014393635030797935899357664859732073 // sqrt(e)*(e - 1)

namespace option {
#pragma warning(disable: 4355)

template<class R = double, class T = double, class S = double, class Z = double>
class gamma : public distribution_base<R,T,S,Z> { // Gamma(a,b)
	typename dual::number_traits<T>::type a_, b_;
public:
	gamma(const T& t, double a = 1, double b = 1)
		: distribution_base(t), a_(GAMMA_ALPHA*a), b_(GAMMA_BETA*b)
	{ }
	~gamma()
	{ }

	// probability density function
	R pdf(const Z& z) const
	{
		return 0; 
	}

	// cumulative distribution function
	R cdf(const Z& z) const
	{
		return exp(lgamma(a_*t_ + 1, b_*exp(z)) - lgamma(a_*t_ + 1) - log(b_)); // = cudf(1, z)
	}

	// kappa(t, s) = log E[exp(s X_t)]
	R cumulant(const S& s) const
	{
		return lgamma(a_*t_ + s) - lgamma(a_*t_) + s*log(b_); // expand around s = 0???
	}

	// cumulant distribution function
	// log E[exp(sX_t) 1(X_t <= x)]
	R cudf(const S& s, const Z& z) const
	{
		return lgamma(a_*t_ + s, b_*exp(z)) - lgamma(a_*t_ + s) - s*log(b_);
	}

	// exp(-kappa(s) + kappa(s, z))
	R N(const S& s, const Z& z) const
	{
		return exp(lgamma(a_*t_ + s, b_*exp(z)) - lgamma(a_*t_ + s));
	}
};

} // namespace option