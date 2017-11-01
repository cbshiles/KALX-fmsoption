// exp_normal_exp.h - analogous to Standard Brownian motion with the t-time normal distribution sqrt(t)*N(0,1) replaced with the sqrt(t)*Normalexp.
#pragma once
#include "option.h"
#include "../fmsdual/normal.h"

namespace option {
#pragma warning(disable: 4355)

template<class R = double, class T = double, class S = double, class Z = double>
class exp_normal_exp : public distribution_base<R,T,S,Z> {
	double a_, alpha_, b_, beta_, den_;
public:
	exp_normal_exp(const T& t, double a, double alpha, double b, double beta)
		: distribution_base(t), a_(a), alpha_(alpha), b_(b), beta_(beta)
	{
		ensure (t > 0);
		ensure (a < 0);
		ensure (alpha > 0);
		ensure (b > 0);
		ensure (beta < 0);
		den_ = normal_pdf(a_)/alpha_ + normal_cdf(b_) - normal_cdf(a_) - normal_pdf(b_)/beta_;
	}

	// cumulative distribution function
	R cdf(const Z& z) const
	{
		Z z_ = z/sqrt(t_);

		if (z_ < a_)
			return exp(alpha_*(z_ - a_))*normal_pdf(a_)/(alpha_*denominator());

		R r = normal_pdf(a_)/alpha_;

		if (z_ < b_)
			return (r + normal_cdf(z_) - normal_cdf(a_))/denominator();

		r += normal_cdf(b_) - normal_cdf(a_);

		return (r + (exp(beta_*(z_ - b_)) - 1)*normal_pdf(b_)/beta_)/denominator();
	}

	// kappa(t, s) = log E[exp(s X_t)]
	// returns log E[exp(s X_t)], where X_t is the normalexp process at time t
	R cumulant(const S& s) const
	{
		// cumulant is integral t infinity
		return 0;//log(normalexp_exp(1999, s*sqrt(t_), a_, alpha_, b_, beta_) ); // ???
	}

	// cumulant distribution function
	// log E[exp(sB_t) 1(B_t <= z)] = s^2t/2 + log P(X_t + st <= z)
	// returns log E[exp(sB_t) 1(B_t <= z)], where B_t is the normalexp process at t
	R cudf(const S& s, const Z& z) const
	{
		return 0;//normalexp_exp(z/sqrt(t_), s*sqrt(t_), a_, alpha_, b_, beta_);
	}

	// exp(-kappa(t,s) + kappa(t,s, z))
	R N(const S& s, const Z& z) const
	{
		return exp(-cumulant(s) + cudf(s, z));
	}
};

} // namespace option