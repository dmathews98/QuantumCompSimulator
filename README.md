# Quantum-Computing-Project

QUANTCOM CIRCUIT SIMULATOR
relevant 18/3/2019

These files can be used to simulate any quantum computing circuit you wish. RunGrover.py is an example of using these files to run grovers algorithm. Remember that n is the number of qubits that have been combined to make a register. The sie of the actual matrix will be 2^n.

To use, import with the line "from dense import *". To use the sparse implementation (not fully working) "from sparse import *". To use the lazy implementation (working even less) "from lazy import *".

In order to create a qubit in a particular state given by an array: Qubit(array)
q.normalise() can be used to normalise if necessary.
In order to create a zero state qubit register of sie 2^n:
Qubit(n)

In order to create Gates:

Some are generalised to act on n sized qubit registers (n=x indicates default size if not specified.
Hadamard: Hadamard(n=1)
Phase shift: Phase(phase shift, n=1)
Identity: Identity(n=1)
Pauli X: Paulix(n=1)
Controlled Not: CNot(n=2)
Controlled Phase: CPhase(phase,n=2)
Swap: Swap(n=2, index to be swapped = 0,index to be swapped = 1)
Controlled versions of any 1 qubit gate can be created: Controlled(gate,n=2)

As a relic, there are some hard coded Combined gates for Grovers:
Diffusion(n)
Oracle(n,fockspace target)

Other gates are not generalised, but can be by using gate.Generalise(n)
1 qubit gates:
V : V()
Pauli Z: PauliZ()
Paulit Y: PauliY()
3 qubit gates:
Toffoli: Toffoli()

To combine gates vertically on a circuit diagram use &
    eg a = Hadamard(x) & PauliZ(y)
    The result will be a combined gate that can act on quantum       registers of size n = x + y

To combine gates horizontally on a circuit diagram use *
    eg a = Hadamard(n) * PauliZ(n)
    The gates must be the same size.
    The result will be a combined gate that enacts both gates when applied to a register of size n

To combine qubits into quantum registers use &
    eg Q = Qubit(x) & Qubit(y)
    The result will be a quantum register of size n = x + y
    Q = Qubit(x + y)

To enact gates on qubits use *
    eg q = Hadamard(n)*Qubit(n)
    the size of the gate and register must be equal.
    The result will be a transformed qubit register
    q = Qubit(n)

While there a fair number of different orders in which gates and qubits can be combined to get the same result that will all work, we recommend the following for speed and clarity:
0. Inspect your circuit diagram, and create the gates you believe will be needed. H = Hadamard(), I = Identity(2), cN = CNot(), PX = PauliX(3)
1. Create a quantum register in 0 state of desired size (count the wires for n): q = Qubit(3)
2. Inspect each column of your circuit diagram and use & to combine each gate that appears top to bottom, using Identity() for blank wires. g1 = H & I, g2 = cN & H, g3 = PX
3. Enact each column on the quantum register in turn:
q = g1*q, q = g2*q, q = g3*q

This implements this circuit:
--- Hadamard --- CNot control --- PauliX
-----------------CNot Target  --- PauliX
-----------------Hadamard     --- PauliX

After applying the gates you want to qubits, the final state of the qubit can be measured using qubit.measure(). This will select a probable final state from the quantum register, in the form of an all 0 array with a 1 in the selected state.

# Using Grover's Algorithm #

Running Grover's has been made mostly user friendly via the terminal. It can only be ran in python3, and done by calling RunGrover.py.

Each function has been described in the class, so there are more details there.

# RunGrover.py #

This file will be set to having check = 'run'. This allows the simulation to run as normal. There are options to run tests by changing check. These are:
- 'test1' which runs for varying register size and graphs the time to run
- 'test2' which runs for varying Fock target and graphs the time to run

# InOut.py #

This controlled the I/O of the programme, with other sections for Shor's etc. Run.Grover.py calls on functions here to determine how to run Grover's, such as size of register, specific or random target value and noise on gates.

# GroverGateWise.py #

This is the main file for governing Grover's algorithm. There are a couple of options to run here.
- In the function 'grover', there are two NOTE sections which describe areas to comment in to run using full preformed Oracle and Diffusion gates, but by default the faster method(sequential application) is used.
- There is a NOTE at the bottom of 'run', which informs of the option to comment in the 'display' function which will display the Oracle and Diffusion matrices used for that simulation. This is done even if sequential application was used as it is the easiest way to see the effects of noise. Using these will increase the runtime.


# shor.py #


This file contains all the functions used to implement Shor's algorithm. To run Shor's algorithm,
import the file or insert a main method, and use shor(N,qubits,noise).

qubits is the number of qubits for the QFT
N should be a non square semi prime that is between 2^(qubits-1) and 2^(qubits)


For testing shors, use either step_test(), noise_test(), or use the results from semi_primes()
