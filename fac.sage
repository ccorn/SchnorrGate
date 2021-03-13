from fpylll import IntegerMatrix, SVP, FPLLL
import sys

def svp(B):
	A = IntegerMatrix.from_matrix(B)
	return SVP.shortest_vector(A)

def is_smooth(x, P):
	y = x
	for p in P:
		while p.divides(y):
			y /= p
	return abs(y) == 1

# Test if a factoring relation was indeed found.
def test_Schnorr(N, n, prec=1000):
	P = Primes()[:n]
	f = list(range(1, n+1))
	shuffle(f)

	# Scale up and round
	def sr(x):
		return round(x * 2^prec)

	N1 = round((N^(1/(n+1))) * (2^prec)) / (2^prec)
	diag = [sr(N1*f[i]) for i in range(n)] + [sr(N1*ln(N))]
	B = diagonal_matrix(diag, sparse=False)
	for i in range(n):
		B[i, n] = sr(N1*ln(P[i]))


	b = svp(B)
	e = [b[i] / diag[i] for i in range(n)]
	en = (b[n] - sum(e[i]*B[i, n] for i in range(n))) / B[n, n]
	assert en in ZZ
	if en == 1:
		print("\nFlipping sign of SVP solution with e_{n+1} == %d" % en)
		e = [-ei for ei in e]
		en = -en	# just for consistency
	elif en != -1:
		print("\nSkipping SVP solution with e_{n+1} == %d" % en)
		return false

	u = 1
	v = 1
	for i in range(n):
		assert e[i] in ZZ
		if e[i] > 0:
			u *= P[i]^e[i]
		if e[i] < 0:
			v *= P[i]^(-e[i])

	return is_smooth(u - v*N, P) 

try:
	bits = int(sys.argv[1])
except:
	bits = 400

try:
	n = int(sys.argv[2])
except:
	n = 47

try:
	trials = int(sys.argv[3])
except:
	trials = 100


FPLLL.set_external_enumerator(None)
print("Testing Schnorr's relation finding algorithm with n=%d on RSA-moduli of %d bits, %d trials"%(n, bits, trials))

successes = 0
for i in range(trials):
	p = random_prime(2^(bits/2), false, 2^(bits/2-1))
	q = random_prime(2^(bits/2), false, 2^(bits/2-1))
	N = p*q
	success = test_Schnorr(N, n)
	successes += success
	print(success, end="\t")
	sys.stdout.flush()

print("\n%d Factoring Relation found out of %d trials"%(successes, trials))
