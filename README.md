# Exceptional cliques code

The code in this repository accompanies the (draft) paper "On exceptional cliques in matrix rings" by Milan Boutros, Ignacio Cascudo, Ronald Cramer, Daniël van Gent and Chaoping Xing, available at https://arxiv.org/abs/2608.24586.

## Files

* **VerifyRecordTable.sage** contains SageMath code that contains the constructions for the exceptional cliques in orders and matrix rings that we list in the paper, and verifies their correctness.

* **NFClique.sage** contains SageMath code that computes the clique number of finite étale algebras.

* **EtaleUpperBound/** is a folder that contains C code that computes the clique number of finite étale algebras, and is optimized for speed.

* **CodingTheoryConstant.sage** contains SageMath code that computes the constant factor in the expression for the asymptotic behaviour of the cliques constructed from coding theory.

Our SageMath code runs on version 10.9.
