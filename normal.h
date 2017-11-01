// normal.h - standard normal distribution
#pragma once

namespace fms {
	template<class I>
	struct normal {
		template<class Z>
		static Z cdf(Z);
	};
} // fms