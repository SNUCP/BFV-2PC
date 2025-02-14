# Ciphertext Simulatable BFV

This repository provides an implementation of "Ciphertext-Simulatable HE from BFV with Randomized Evaluation"(https://eprint.iacr.org/2025/203), based on HIENAA library (unreleased).

Before you run the code, please make sure to install the following packages : ChaChaCiphers, Primes, Nemo, NormalForms, LinearAlgebra, DoubleFloats.
To install them, you can open the REPL and type the following commands.

<pre>
<code>
]
add ChaChaCiphers
add Primes
add Nemo
add NormalForms
add LinearAlgebra
add DoubleFloats
</code>
</pre>

To run the test code, type the following command in the terminal.

<pre>
<code>
julia distance.jl
</code>
</pre>
