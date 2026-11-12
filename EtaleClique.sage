# Requires sage 10.4
import itertools
import sys

### Graphics
def print_progress_bar(index,total,score):
    n_bar = 20 
    progress = index / total
    sys.stdout.write('\r')
    sys.stdout.write(f"[{'=' * int(n_bar * progress):{n_bar}s}] {index}/{total}: {score}")
    sys.stdout.flush()

### Enumerate etale algebras
def partitions(n):
    if n <= 0:
        return []
    def partitions_aux(target,least):
        if target == 0:
            return [[]]
        P = []
        for i in range(least, target+1):
            for p in partitions_aux(target-i, i):
                P.append([i] + p)
        return P
    return partitions_aux(n, 1)

def rings_to_check(n,p):
	for P in partitions(n):
		yield P, itertools.product(*(GF(p^k) for k in P))

### Auxiliary functions
def are_connected(x,y):
	det = 1
	for xi,yi in zip(x,y):
		det *= (xi-yi).norm()
	return det in [-1,1]

# if c is None, returns the maximal clique number among all etale F_p-algebras of dimension n
# if c is a positive integer, returns whether a clique of size c exists among these algebras
def check_rings(n,p,c=None,verbose=False):
	best_clique_number = 0
	# enumerate all rings R
	for s, R in rings_to_check(n,p):
		zero = [GF(p^x)(0) for x in s]
		one = [GF(p^x)(1) for x in s]
		if verbose:
			print("\r::: R =", ' x '.join("GF(" + str(p) + "^" + str(x) + ")" for x in s), "             ")
		# compute all exceptional cliques {0,1,x} in R 
		U = [ x for x in R if are_connected(x,zero) and are_connected(x,one) ]
		if 2+len(U) < c:
			if verbose:
				print(": too few exceptional units")
			continue
		if verbose:
			print(": # exceptional units =", len(U))
		# compute the equivalence classes under Aut(R) and the actions of x -> 1/x and x -> 1-x
		E = { frozenset(itertools.chain((xi.minpoly() for xi in x), ((1/xi).minpoly() for xi in x), ((1-xi).minpoly() for xi in x), ((1/(1-xi)).minpoly() for xi in x), ((xi/(xi-1)).minpoly() for xi in x), (((xi-1)/xi).minpoly() for xi in x))) : x for x in U }
		if verbose:
			print(": # equivalence classes =", len(E))
		# choose some representative x
		local_best_clique_number = 0
		for i, (xc, x) in enumerate(E.items()):
			# compute all exceptional cliques {0,1,x,y}
			L = [y for y in U if are_connected(x,y)]
			if verbose:
				print_progress_bar(i+1, len(E), str(local_best_clique_number)+" "+str(len(L)))
			G = Graph({ y : [ z for z in L if are_connected(y,z) ] for y in L })
			clique_number = G.clique_number() + 3
			if c != None:
				if clique_number >= c:
					if verbose:
						print("\n::: a clique of size", c, "exists")
					return True
			else:
				if clique_number > local_best_clique_number:
					local_best_clique_number = clique_number
			U = [ y for y in U if (yi.minpoly() for yi in y) not in xc ]
		if verbose and c == None:
			print("\n: largest clique =",local_best_clique_number)
		best_clique_number = max(best_clique_number,local_best_clique_number)
	if c != None:
		if verbose:
			print("::: no cliques of size", c, "exist")
		return False
	else:
		return best_clique_number
