# The numeric constant factor for the nlog(n) clique family
def compute_constant(d,q):
	epsilon = 2 if q % 2 == 0 else 1
	delta = (q/d)/RR(log(q))
	for i in range(20):
		delta /= (1+epsilon/q^(2^i))
	return delta
