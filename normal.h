// normal.h - standard normal distribution
#pragma once
#include <cmath>

namespace fms {
	template<class I>
	struct normal {

		// cumulative distribution
		template<class Z>
		static Z cdf(Z);

		// inverse cdf
		template<class Z>
		static Z inv(Z);

		// probability density
		template<class Z>
		static Z pdf(Z z)
		{
			return exp(-z*z/2)/M_SQRT_2PI;
		}
	};
} // fms