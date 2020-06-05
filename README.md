# Shapes of number fields

This repository contains one main file: shape.sage. This is a sage script file defining functions for computing the shape, the trace-zero form, and the unit shape of a number field.

The shape of a number field *K* is an invariant that was introduced in David Terr's PhD thesis. It basically measures the 'shape' of the ring of integers of *K* as a lattice. The trace-zero form is a related, but different invariant. Note that some authors use the term 'shape' to refer to the trace-zero form (e.g. Bhargava–Shnidman), and that some authors take a slightly different definition of these two invariants (e.g. Mantilla-Soler). The definitions in use here follow Terr's original definition. See e.g. Section 2 of my article *Equidistribution of shapes of complex cubic fields of fixed quadratic resolvent* ([arXiv:1907.07209 [math.NT]](https://arxiv.org/abs/1907.07209)) for more information.

The unit shape of a number field *K* is the 'shape' of its unit lattice in multiplicative Minkowski space. Very little is known about this invariant. See e.g. [arXiv:1910.12105 [math.NT]](https://arxiv.org/pdf/1910.12105.pdf) for a recent result.