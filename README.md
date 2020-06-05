# Shapes of number fields

This repository contains one main file: shape.sage. This is a sage script file defining functions for computing the shape, the trace-zero form, and the unit shape of a number field.

The shape of a number field *K* is an invariant that was introduced in David Terr's PhD thesis. It basically measures the 'shape' of the ring of integers of *K* as a lattice. The trace-zero form is a related, but different invariant. Note that some authors use the term 'shape' to refer to the trace-zero form (e.g. Bhargava–Shnidman), and that some authors take a slightly different definition of these two invariants (e.g. Mantilla-Soler). The definitions in use here follow Terr's original definition. See e.g. Section 2 of my article *Equidistribution of shapes of complex cubic fields of fixed quadratic resolvent* ([arXiv:1907.07209 [math.NT]](https://arxiv.org/abs/1907.07209)) for more information.

The unit shape of a number field *K* is the 'shape' of its unit lattice in multiplicative Minkowski space. Very little is known about this invariant. See e.g. [arXiv:1910.12105 [math.NT]](https://arxiv.org/pdf/1910.12105.pdf) for a recent result.

## Sample usage

### Shapes
The shape of the real subfield of the 7th cyclotomic field. This field is totally real and so by default its shape is returned as a matrix with coefficients in the integers (the function `shape_of_a_number_field(K)` detects if `K` is totally real). (All shapes of Galois cubic fields are hexagonal lattices).

```
sage: K = NumberField(x^3-x^2-2*x+1, 'a')
sage: shape_of_a_number_field(K)
[42 21]
[21 42]
```

Setting `check` to `False` skips checking if `K` is totally real and defaults to doing computations in the algebraic reals and returns a matrix over that field.
```
sage: shape_of_a_number_field(K, check=False)
[42.00000000000000? 21.00000000000000?]
[21.00000000000000? 42.00000000000000?]
```

If you know that your field will be totally real, you can skip the check and use the following:
```
sage: K = NumberField(x^3 - x^2 - 2*x + 1, 'a')
sage: shape_of_a_totally_real_number_field(K)
[42 21]
[21 42]
```

The shape of a non-totally real number field will always be returned as a matrix over the algebraic reals:
```
sage: K = NumberField(x^3 - 10, 'a')
sage: shape_of_a_number_field(K)
[ 78.55780720179485? -36.78350769927984?]
[-36.78350769927984? 120.33210670430986?]
```

### Trace-zero forms
The trace-zero form of a totally real number field is the same as its shape. For a non-totally real number field, it will *always* be different as it is no longer a positive-definite quadratic form. It does however always have integer coefficients:
```
sage: K = NumberField(x^3 - 10, 'a')
sage: trace_zero_form(K)
[  60   30]
[  30 -120]
```

### Unit shapes
The unit shapes involve logarithms of algebraic numbers and are therefore returned with entries in the real numbers. Here is the unit shape of the real subfield of the 7th cyclotomic field (all unit shapes of Galois cubic fields are hexagonal lattices).
```
sage: K.<a> = NumberField(x^3 - x^2 - 2*x + 1)
sage: unit_shape(K)
[ 1.05090936424514477667765209089664901908211009202358957014 0.525454682122572388338826045448324509541055046011794785069]
[0.525454682122572388338826045448324509541055046011794785069  1.05090936424514477667765209089664901908211009202358957014]
```

Here is an example of the biquadratic field **Q**(&radic;2, &radic;3) whose unit shape is a right rectangular prism:
```
sage: K = NumberField(x^4 - 10*x^2 + 1, 'a')
sage: unit_shape(K)
[     1.73437810227263615040825812545843864708057890217645004209  3.47294036626874596392295926713361173618860793031990895364e-56 -1.27447352890596182162310431821416944447288364415409502886e-57]
[ 3.47294036626874596392295926713361173618860793031990895364e-56      3.10727759958278392600229884035767561842823125723981963686 -3.21804566048755359959833840349077784729403120148908994787e-56]
[-1.27447352890596182162310431821416944447288364415409502886e-57 -3.21804566048755359959833840349077784729403120148908994787e-56      5.25524295960704856821635998862054522346506619464633121785]
```
