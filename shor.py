import numpy as np
from dense import *
#from sparse import *
import InOut as IO
import time
import matplotlib.pyplot as plt
from fractions import Fraction


"""Code to implement Shor's algorithm for polynomial time factoring of numbers.
Broken up into various subroutines (some quantum some classical)
"""

#----------------------------------------------------
#--- Quantum subroutines --------------- AR ---------
#----------------------------------------------------

#--- Helper functions for QFT ---

def Flip(n,noise=0):
	"""Applies multiple swap gates to flip the order of all qubits. Optional noise
	"""
	if noise==0:
		g = Identity(n)
		for i in range(n//2):
			if i!=(n-1-i):
				g = g*Swap(n,i,n-1-i)
		return g
	else:
		g = Identity(n)
		for i in range(n//2):
			if i!=(n-1-i):
				g = g*Noisy(Swap(n,i,n-1-i),noise)
		return g


def R(n,noise=0):
	"""Function that calls controlled phase gate, with optional noise
	"""
	if noise==0:
		return CPhase(-np.pi*2/(2.0**n))
	else:
		return Noisy(CPhase(-np.pi*2/(2.0**n)),noise)

def H(n=1,noise=0):
	"""Function that calls Hadamard gate of size n, with optional noise
	"""
	if noise==0:
		return Hadamard(n)
	else:
		return Noisy(Hadamard(n),noise)
def I(n=1,noise=0):
	"""Function that calls Identity gate of size n, with optional noise
	"""
	if noise==0:
		return Identity(n)
	else:
		return Noisy(Identity(n),noise)

def S(n,a,b,noise=0):
	"""Function that calls Swap gate of size n, with optional noise
	"""
	if noise==0:
		return Swap(n,a,b)
	else:
		return Noisy(Swap(n,a,b),noise)

#--- QFT defined recursively ---


def QFT(N,noise=[0,0,0,0]):
	"""Function to combine gates to form a Quantum Fourier Transform (QFT) on N qubits.
	Defined recursively. Returns Identity gate if 4 QFTs are multiplied in sequence
	noise vector applies noise to Swap, Phase, Hadamard and Identity gates
	"""
	def _QFT(n):	
		if n==2:
			#Base case
			return ((H(noise=noise[2])&I(noise=noise[3]))*(R(2,noise=noise[1]))*(I(noise=noise[3])&H(noise=noise[2])))#*Swap(2,0,1))
		else:
			g1 = _QFT(n-1)&I(noise=noise[3])
			g2 = I(n,noise=noise[3])
			for i in range(n-2):
				g2 = g2*S(n,i,n-2,noise=noise[0])*(I(n-2,noise=noise[3])&R(n-i,noise=noise[1]))*S(n,i,n-2,noise=noise[0])
			g3 = (I(n-2,noise=noise[3])&R(2,noise=noise[1]))*(I(n-1,noise=noise[3])&H(noise=noise[2]))
			return (g1*g2*g3)
	#Some ambiguity whether Flip gate is needed here. Appears to work with or without
	return _QFT(N)*Flip(N,noise=noise[0])







def iQFT(n,noise=[0,0,0,0]):
	"""
	Inverse QFT, defined using the maths that 
	QFT^4 = Identity => (QFT^3)*QFT = Identity => inverse of QFT = QFT^3
	"""
	f = QFT(n,noise)
	return f*f*f




#-------------------------------------------------------------------------
#--- Classical subroutines -----------------------------------------------
#-------------------------------------------------------------------------



def GCD(x,y):
	"""Euclidean algorithm for greatest common divisor between 2 numbers
	"""
	if x==y:
		return x
	while(y): 
		x,y = y,x%y 
	return x

def extendedGCD(x,y):
	"""Extended euclidean algorithm that returns array of remainders
	that arise between 2 numbers during the process of finding the GCD.
	"""
	fracs = []
	if x==y:
		return fracs
	while(y):
		fracs.append(x//y)
		x,y = y,x%y
	return np.array(fracs)

def isPrime(n):
	"""Method to check if a given number is prime. Used to check that Shors isn't
	given a prime number to factorise
	"""
	if n <= 3:
		return n>1
	elif (n%2==0) or (n%3==0):
		return False
	else:
		i = 5
		while i*i <= n:
			if (n%i==0) or (n%(i+2)==0):
				return False
			i = i+6
		return True


def modexp(b,N,qubits):
	"""Evaluates the modular exponential function (base b, mod N) across an array.
	"""
	def _modexp(b,e,N):
		t = 1
		for i in range(e):
			t = b*t
		return t%N
	
	xs = np.arange(0,2**qubits,1)
	g = np.vectorize(lambda x:_modexp(b,x,N))
	f = lambda x:(g(x))
	return f(xs)/np.sqrt(np.sum(np.square(xs).astype(float)))




def shor(N,qubits,noise=[0,0,0,0],purely_quantum=False,verbose=False):
	"""Shors algorithm for number factorisation. N is non square semiprime number to factorise.
	qubits is number of qubits to use for QFT. purely_quantum defines whether to allow
	initial guess being correct to terminate run, used mainly for checking the rest of the
	algorithm works, as often the initial guess being correct is more likely than the QFT
	factorising (especially for smaller numbers). verbose triggers various print statements

	"""
	
	assert N%2!=0, "N must be odd"
	assert isPrime(N)==False, "Cannot factorise prime number"+str(N)

	guess = np.random.randint(2,(N-1))
	p = 1
	counter = 0
	found = False

	qubit_profile = np.zeros(2**qubits)

	
	while found==False:#p%2==1 or (guess**(p/2))%N != (N-1):
		counter = counter +1
		guess = np.random.randint(2,(N-1))
		
		divisor = GCD(N,guess)
		if verbose:
			print("Does "+str(guess)+" factor "+str(N)+" ?")

		if divisor!=1:
			
			#If gcd(N,guess) is not 1, then guess is a factor of N. Lucky
			if verbose:
				print("Divisor "+str(divisor)+" factorises "+str(N)+", found by guessing")
			#found = True
			if purely_quantum:
				pass
			else:
				return [divisor,divisor,counter]

			#pass

		#Run modular exponentiation with base = guess, mod = N
		qreg = Qubit(modexp(guess,N,qubits))
		#print(np.sum(qreg.ret_mod()))
		#Apply QFT to qubit register
		qreg = iQFT(qubits,noise)*qreg
		
		#qubit_profile = np.vstack((qubit_profile,qreg.ret_mod()))
		#Measure qubit register to estimate frequency

		#IO.Hist(qreg)
		qreg.measure_cheat()
		freq = int(str(qreg.split_register()),2)
		
		#Transform frequency to nearest discrete value
		sample_freq = Fraction(freq,2**qubits)
		p = sample_freq.limit_denominator(N).denominator
		if verbose:
			print("Period guess "+str(p))

		#if period is odd or is trivial square root, start again

		if p%2==0 and ((guess**p/2)%N!=(N-1) and (guess**p/2)%N!=1):
			g1 = GCD(guess**p/2+1,N)
			g2 = GCD(guess**p/2-1,N)
			#print(g1)
			#print(g2)

			if g1!=1:
				if verbose:
					print(str(g1)+" factorises "+str(N)+", found with QFT")
				found = True
				#return
			if g2!=1:
				if verbose:
					print(str(g2)+" factorises "+str(N)+", found with QFT")
				found = True
			if g1!=1 or g2!=1:
				#print(str(counter)+" steps")
				#return counter
				#plt.matshow(qubit_profile)
				#plt.show()
				return [g1,g2,counter]
	return shor(N,qubits)
	

	



#---------------------------------------------------------
#--- Methods for testing shors ---------------------------
#---------------------------------------------------------

def semi_primes(bits,bits_lower = 1):
	"""Randomly generates a product of 2 primes, n, such that
	2**bits_lower < n < 2**bits
	Also returns the 2 prime factors. Very useful for testing shors
	"""
	#print(list(filter(isPrime,range(2**bits))))
	primes = np.array(list(filter(isPrime,range(2**bits))))[1:]
	prime_products = np.outer(primes,primes)

	res = 0
	while (res >2**bits) or (res < 2**(bits_lower)):
		x = np.random.randint(0,len(primes))
		if x==1:
			y=0
		elif x==0:
			y=1
		else:
			y = np.random.randint(0,x-1)
		res = prime_products[x,y]
	return res,primes[x],primes[y]


def step_test(lower,upper,its=10,noise=[0,0,0,0]):
	"""Method for testing shors on randomly generated prime products.
	Returns mean number of QFTs applied and runtimes, with variances. 
	upper and lower give upper and lower bounds to number of qubits in QFT. 
	its defines how many numbers to attempt to factorise at each size
	"""
	hits = 0
	misses = 0
	
	steps = np.zeros((its,len(range(lower,upper))))
	times = np.zeros((its,len(range(lower,upper))))
	for bits in (range(lower,upper)):
		print("Running tests for "+str(bits)+" Qubits")
		for x in range(its):
			print(x)
			coprime,prime1,prime2 = semi_primes(bits+1,bits-2)
			t1 = time.time()
			res = shor(coprime,bits,noise)
			t2 = time.time()
			steps[x,bits-lower] = res[2]
			times[x,bits-lower] = t2-t1
			if (res[0]==prime1) or (res[0]==prime2) or (res[1]==prime1) or (res[1]==prime2):
				hits = hits+1
			else:
				misses = misses +1
	ms = np.mean(steps,axis=0)
	ers = np.std(steps,axis=0)

	m_time = np.mean(times,axis=0)
	ers_time = np.std(times,axis=0)

	print(str(hits)+" succesful factorings")
	print(str(misses)+" failures")

	return ms,ers,m_time,ers_time


def noise_tests(min_qubits,max_qubits,its):

	max_qubits = max_qubits+1
	noise_hadamard=[0,0,0.3,0]
	noise_phase = [0,0.3,0,0]
	ms,errs,ts,t_errs = step_test(min_qubits,max_qubits,its)
	hms,herrs,hts,ht_errs = step_test(min_qubits,max_qubits,its,noise_hadamard)
	pms,perrs,pts,pt_errs = step_test(min_qubits,max_qubits,its,noise_phase)
	xs = range(min_qubits,max_qubits)
	plt.errorbar(x=xs,y=ms,yerr=errs,fmt='o',label="No noise",capsize=10)
	plt.errorbar(x=xs,y=hms,yerr=herrs,fmt='o',label="Hadamard noise",capsize=10)
	plt.errorbar(x=xs,y=pms,yerr=perrs,fmt='o',label="Phase noise",capsize=10)
	plt.legend()
	plt.xlabel("Number of Qubits")
	plt.ylabel("Number of QFT applications")
	plt.show()



	plt.errorbar(x=xs,y=ts,yerr=t_errs,fmt='o',label="No noise",capsize=10)
	plt.errorbar(x=xs,y=hts,yerr=ht_errs,fmt='o',label="Hadamard noise",capsize=10)
	plt.errorbar(x=xs,y=pts,yerr=pt_errs,fmt='o',label="Phase noise",capsize=10)
	plt.legend()
	plt.xlabel("Number of Qubits")
	plt.ylabel("Time taken (s)")
	plt.show()	




def main():




	#noise_tests(3,6,20)
	noise_tests(3,6,5)
	#print(semi_primes(7,4))
	#print(shor(7*5,5,noise,verbose=True))
	
	#print(isPrime(209))
	#for x in range(10):
	#	print(semi_primes(8))
	#ls = filter(isPrime,range(1000))
	#print(ls)
	#a = Fraction(1.234567).limit_denominator(1000)

	#print(a)

	#print(extendedGCD(416,93))
	#print(continued_fraction(1,93,416))
	#q1 = Qubit(8)
	#print(q1.ret_mod())
	#q = (Hadamard(3)&Identity(5))*q1
	#qreg = Qubit(modexp(41,37,6))
	#IO.Hist(qreg)
	#noise = [0,0,0,0]
	#ift = iQFT(3,noise)
	#ft = QFT(3,noise)
	#ift = ft*ft*ft
	#IO.Display(QFT(4,noise)*QFT(4,noise)*iQFT(4,noise)*iQFT(4,noise))
	#IO.Display(ft)
	#IO.Display(ift)
	#IO.Display(ft*ift)
	#IO.Hist(ft*ift*qreg)
	
	#IO.Display(ft)
	#IO.Hist(ft*qreg)
	
	"""
	qreg = Qubit(modexp(13,29,6))
	IO.Hist(qreg)
	noise = [0,0,0,0]
	ft = QFT(6,noise)
	IO.Display(ft)
	IO.Hist(ft*qreg)
	"""
	

main()
