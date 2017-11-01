// implied.h - implied volatilty
// Copyright (c) 2011 KALX, LLC. All rights reserved. No warranty is made.
#pragma once
#include <functional>
#include "../fmsdual/number.h"
#include "../fmsroot/newton.h"

namespace option {

	// |f(s) - p| < eps with initial guess s.
	template<class T>
	inline T 
	root1d(std::function<dual::number<T,3>(dual::number<T,3>)> F, const T& p, T s, const T& eps)
	{
		dual::number<T,3> f;

		
		for (f = F(dual::number<T,3>(s, 1)) - p; fabs(f[0]) > eps; f = F(dual::number<T,3>(s, 1)) - p)
			s = root::_1d::newton_2(s, f.data());

		return s;
	}

} // namespace option