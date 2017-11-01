Option pricing for arbitrary distributions

This project makes it simple to price European puts, calls, and binary options. In addition, numerically stable calculations of greeks of all orders can be computed using the fmsdual project.

If X is a random variable then the cumulant of X is κ(s) = log E[exp(s X)]. If (Xt)t≥0 is a stochastic process the cumulant of Xt is κt(s) = log E[exp(s Xt)]. If (Xt) is a Lévy process then κt(s) = t κ(s), but that is not used here.

Define the cumulant distribution function by κt(s, z) = log E[exp(s Xt) 1(Xt ≤ z)].

The Fischer Black model gives the foward price of a put option as E[max{k - F, 0}], where k is the strike and F is the risk-neutral value of the underlying at expiration. If F = f exp(-κt(s) + s Xt), then E[F] = f. The put price is

E[max{k - F, 0}]	=	k E[1(F ≤ k)] - E[F 1(F ≤ k)]
 	=	k P(F ≤ k) - f E[exp(-κt(s) + sXt) 1(F ≤ k)]
 	=	k P(Xt ≤ z) - f P*(Xt ≤ z),
where z = (log k/f + κt(s))/s and P* is the Esscher transformation of P, i.e., dP*/dP = exp(sXt)/E[exp(sXt)] = exp(-κt(s) + sXt).

If Xt is standard Brownian motion, then κ(t, s) = s2t/2 and κt(s, z) = κt(s) - log N((z - st)/√t), where N is the standard normal cummulative distribution. The above reduces to

k N(-d2) - f N(-d1),
for the standard values of d1 and d2.