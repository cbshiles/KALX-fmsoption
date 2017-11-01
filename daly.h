// daly.h
// From: bill.daly@tradition-ny.com (Bill Daly)
// Message-ID: <4fa0d5e.0402241117.53c8da2a@posting.google.com>
// Message-ID: <4fa0d5e.0402181449.6dd203db@posting.google.com>
// standard normal cumulative distribution function
#pragma once
#include <cmath>
#include "normal.h"

namespace fms {
	class daly {};
	template<>
	struct normal<daly> {
		template<class Z>
		static Z cdf(Z x)
		{
			Z s=x,t=0,b=x,x2=x*x,i=1;
	
			while(s!=t)
				s=(t=s)+(b*=x2/(i+=2));

			return  static_cast<Z>(.5+s*exp(-.5*x2-.91893853320467274178L));
		}
	};
} // fms