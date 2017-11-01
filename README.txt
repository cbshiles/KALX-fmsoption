GENERAL

Option pricing for arbitrary distributions using modern C++.
Uses dual numbers to calculate greeks. http://fmsdual.codeplex.com

Basic idea:

	put p(k,t); // pays max{k - F(t), 0} at t
	model bm(brownian(f,s)); // F(t) = f exp(-s^2 t/2 + s B(t))

	double v = value(p, bm); // value of put using Brownian motion	

FINANCE/MATHEMATICS

For any random variable F we have

 E[max{k - F, 0}] = k P(F <= k) - E[F] P_(F <= k),
  
where dP_/dP = F/E[F] by a trivial application of the Esscher transform.

Let (X(t)) be any process (but think Brownian motion).
Define F(t) = f exp(-kappa(s, t) + s X(t)) where kappa(s,t) = log E[s X(t)].
(For Brownian motion kappa(s,t) = s^2 t/2.)
(For Levy processes kappa(s,t) = t kappa(s,1).)

Note F(t) <= k iff X(t) <= (kappa(s,t) + log k/f)/s =: z(f,s,t).

SOFTWARE
Need to know auto, decltype, and trailing return types.
