// option.h - Black (forward) option valuation for arbitrary distributions.
// Copyright (c) 2011-2013 KALX, LLC. All rights reserved. No warranty is made.
#pragma once
#ifndef ensure
#include <cassert>
#define ensure(x) assert(x)
#endif
#include <algorithm>

namespace fms {
namespace option {

	// X implements kappa, cdf, and cdf_
	template<class F, class S, class K, class T, class X>
	inline auto put(F f, S s, K k, T t, X x) -> decltype(f+s+k+t)
	{
		typedef decltype(f+s+k+t) Z;

		ensure (f >= 0);
		ensure (s >= 0);
		ensure (k >= 0);
		ensure (t >= 0);

		// special cases at machine epsilon
		if (1 + f == f)
			return k;
		if (1 + s == s || 1 + t == t)
			return std::max<F>(k - f, 0.);
		if (1 + k == k)
			return Z(0.);

		Z z = (x.kappa(s,t) + log(k/f))/s;

		return k*x.cdf(z,t) - f*x.cdf_(z,t,s);
	}

	template<class F, class S, class K, class T, class X>
	inline auto call(F f, S s, K k, T t, X x) -> decltype(f+s+k+t)
	{
		return f - k + put(f, s, k, t, x); // put call parity
	}

	template<class F, class S, class K, class T, class X>
	inline auto value(F f, S s, K k, T t, X x) -> decltype(f+s+k+t)
	{
		return k < 0 ? put(f, s, -k, t, x) : call(f, s, k, t, x);
	}
	/*
	template<class K, class T>
	struct put {
		K k;
		T t;
	};
	template<class K, class T>
	struct call {
		K k;
		T t;
	};

	template<class F, class S, class X>
	struct model : public X {};

	template<class K, class T, template<class,class> O<K,T>, class F, class S, class X, template<class,class,class> M<F,S,X>>
	inline T value(const O<K,T>& o, const M<F,S,X>& m);
	*/
} // option
} // fms
