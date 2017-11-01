// brownian.h - Brownian motion.
// Copyright (c) 2013 KALX, LLC. All rights reserved. No warranty is made.
#pragma once
#include "normal.h"

namespace fms {

	// normal<I>::cdf is any standard normal distribution
	template<class I, class Z = double, class S = double, class T = double>
	struct brownian {
		static auto kappa(S s, T t) -> decltype(s*t) 
		{
			return s*s*t/2;
		}
		static Z cdf(Z z, T t)
		{
			return normal<I>::cdf<Z>(z/sqrt(t));
		}
		static Z cdf_(Z z, T t, S s)
		{
			return cdf(z - s*t, t);
		}
	};

} // fms