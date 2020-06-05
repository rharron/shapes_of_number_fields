from sage.rings.number_field.number_field_base import is_NumberField
from sage.libs.pari.convert_sage import gen_to_sage

def shape_of_a_number_field(K, check = True):
    r"""
    Return the shape of the number field as a Gram matrix.

    The shape of a number field of degree d is by definition the equivalence
    class of the (d-1)-dimensional lattice given by the orthogonal projection of
    the ring of integers of K onto the trace zero space (i.e. onto the
    orthogonal complement of 1) (equivalence means up to change of basis,
    orthogonal transformation, and homothety). The output of this function is
    the Gram matrix of the lattice. By default, this function checks if K is
    totally real in which case it returns the Gram matrix over the integers.
    Otherwise (or if ``check`` is set to False), the algebraic reals are used.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3-x^2-2*x+1)
        sage: shape_of_a_number_field(K)
        [42 21]
        [21 42]
        sage: shape_of_a_number_field(K, check=False)
        [42.00000000000000? 21.00000000000000?]
        [21.00000000000000? 42.00000000000000?]

        sage: K.<a> = NumberField(x^4 + 4*x^2 + 2)
        sage: shape_of_a_number_field(K)
        [128.0000000000000?                  0            0.?e-16]
        [                 0 128.0000000000000?                  0]
        [           0.?e-16                  0 128.0000000000000?]
    """

    if K.absolute_degree() == 1:
        return matrix(ZZ, 0)

    if K.absolute_degree() == 2:
        return Matrix(ZZ, 1, [1])

    if check:
        if K.is_totally_real():
            return shape_of_a_totally_real_number_field(K)

    #Otherwise do it with AA and QQbar
    if K.is_relative():
        KK = K.absolute_field(K.variable_names()[0] * 2)
    else:
        KK = K
    B = [KK(aa) for aa in KK.pari_zk()]
    Bperp = []
    for a in range(1, KK.absolute_degree()):
        Bperp.append(_perp(B[a]))
    sigmas = list(KK.embeddings(AA))
    iotas = list(KK.embeddings(QQbar))
    taus = []
    while len(iotas) > 0:
        iota = iotas.pop()
        aa = iota.im_gens()[0]
        if aa.imag() == 0:
            continue
        taus.append(iota)
        aaconj = aa.conjugate()
        for iota2 in iotas:
            if iota2.im_gens()[0] == aaconj:
                iotas.remove(iota2)
                break
    B = [minkowski_vector(b, sigmas, taus) for b in Bperp]
    IP = diagonal_matrix([1]*len(sigmas) + [2] * (2*len(taus)))
    Gram = Matrix(AA, KK.degree() - 1, [B[i] * IP * B[j] for i in range(len(B)) for j in range(len(B))])
    return Gram

def shape_of_a_totally_real_number_field(K, check = False):
    """
    Return the shape of the totally real number field K.

    The shape of a totally real number field is a quadratic form over the
    rational numbers and the output of this function is the Gram matrix (over
    QQ). The computations are purely algebraic as the Minkowski inner product is
    simply the trace form. By default, this function does not check that the
    given field is totally real. To force a check, set "check" to True.

    """

    if check:
        if not K.is_totally_real():
            raise ValueError("Field passed to shape_of_a_totally_real_number_field is not totally real.")
    return trace_zero_form(K)

def trace_zero_form(K):
    r"""
    Return the Gram matrix of the trace zero form of the number field ``K``.

    INPUT:

    - ``K`` -- a number field or an irreducible polynomial over the rationals

    OUTPUT:

    - The Gram matrix of the trace form projected to the trace zero space
    with the definition that always gives an integral quadratic form.
    """
    if is_NumberField(K):
        T2 = K.pari_nf()[4][3]
    else: #assumes K is a polynomial over QQ or anything that can be converted to pari and passed to nfinit
        K = pari(K)
        T2 = K.nfinit()[4][3]
    G = gen_to_sage(T2)
    return _gram_to_perp_gram(G)

def unit_shape(K, algorithm='pari'):
    r"""
    Return the unit shape of the number field K as a Gram matrix.

    The unit shape of a number field is the equivalence class of its unit
    lattice under rotations, reflections, and scaling. Its unit lattice is the
    lattice obtained by embedding its group of units into CC^r using
    multiplicative Minkowski theory (where r is the rank of its units). See
    Section I.5 and I.7 of Neukirch's Algebraic Number Theory for more details.

    INPUT:

    - ``K`` -- a number field
    - ``algorithm`` -- (Default: 'pari') This can be either 'pari' or 'sage'.
    The 'pari' value is much quicker as it uses data precomputed by pari when
    ``K`` is created.

    OUTPUT:

    - The Gram matrix of the log embeddings of a set of fundamental units of ``K``.

    EXAMPLES::

        sage: K.<a> = NumberField(x^3-x^2-2*x+1)
        sage: unit_shape(K)  #The Gram matrix of a hexagonal lattice
        [ 1.05090936424514477667765209089664901908211009202358957014 0.525454682122572388338826045448324509541055046011794785069]
        [0.525454682122572388338826045448324509541055046011794785069  1.05090936424514477667765209089664901908211009202358957014]

        sage: K.<a> = NumberField(x^5-2)
        sage: G = unit_shape(K); G
        [ 15.0129023175321854981103129410192456961616713603540086007 -4.30063531991780989639517637889332226987381497518212054803]
        [-4.30063531991780989639517637889332226987381497518212054803  5.90327337378657885087428238031757753150997155164921699972]
        sage: r1, r2 = K.signature()
        sage: (K.regulator() * (r1 + r2).sqrt()).n(digits=10)
        8.374353847
        sage: G.det().sqrt().n(digits=10)
        8.374353847
    """
    if algorithm == 'pari':
        return unit_shape_pari(K)
    elif algorith == 'sage':
        return unit_shape_sage(K)
    raise ValueError("Parameter algorithm in unit_shape must be either \'pari\' or \'sage\'.")

def _gram_to_perp_gram(G):
    """
    Input an nxn Gram matrix where the first row and column are assumed to
    correspond to the number 1 and output an (n-1)x(n-1) matrix that is the
    Gram matrix of the perps of the remaining elements.
    """
    n = G[0, 0]
    Gperp = G.submatrix(1, 1)
    for i in range(n-1):
        for j in range(i, n-1):
            Gperp[i,j] = n^2 * Gperp[i, j] - n * G[0, i+1] * G[0, j+1]
            if j != i:
                Gperp[j,i] = Gperp[i,j]
    return Gperp

def minkowski_vector(alpha, sigmas, taus):
    r"""
    INPUT:
        -``alpha`` -- an element of some numebr field K
        -``sigmas`` -- a list (or tuple) of the embeddings of K into AA
        -``taus`` -- a list (or tuple) of the non-real places of K (as embeddings into QQbar)

    OUTPUT:
        -A vector over AA that is the image of ``alpha`` under the Minkowski embedding of K
    """
    n = len(sigmas) + 2 * len(taus)
    V = AA^n
    return V([sigma(alpha) for sigma in sigmas] + sum([[tau(alpha).real(), tau(alpha).imag()] for tau in taus], []))

def multiplicative_minkowski_vector(u, sigmas, taus):
    r"""
    INPUT:
        -``alpha`` -- an element of some numebr field K
        -``sigmas`` -- a list (or tuple) of the embeddings of K into AA
        -``taus`` -- a list (or tuple) of the non-real places of K (as embeddings into QQbar)

    OUTPUT:
        -A vector over RR that is the image of ``alpha`` under the multiplicative Minkowski
        embedding of K
    """
    m = len(sigmas) + len(taus)
    v = [sigma(u).abs().log() for sigma in sigmas] + [(tau(u).abs()^2).log() for tau in taus]
    return (RR^m)(v)

def unit_shape_sage(K):
    r"""
    Return the Gram matrix of the unit lattice of ``K`` using sage's K.units().
    """
    us = K.units(proof=False)
    sigmas = list(K.embeddings(RR))
    iotas = list(K.embeddings(CC))
    r = len(us)
    taus = []
    while len(iotas) > 0:
        iota = iotas.pop()
        aa = iota.im_gens()[0]
        if aa.imag() == 0:
            continue
        taus.append(iota)
        aaconj = aa.conjugate()
        for iota2 in iotas:
            if iota2.im_gens()[0] == aaconj:
                iotas.remove(iota2)
                break
    B = [multiplicative_minkowski_vector(u, sigmas, taus) for u in us]
    return Matrix(RR, r, [[u.dot_product(v) for v in B] for u in B])

def unit_shape_pari(K):
    r"""
    Return the Gram matrix of the unit lattice of the number field ``K`` using
    the data are stored in the underlying pari bnf structure of ``K``.
    """
    r1, r2 = K.signature()
    r = r1 + r2 -1
    bnf = K.pari_bnf(units=False)
    UL = [[bnf[2][j][i].sage().real() for i in range(r+1)] for j in range(r)]
    RR = UL[0][0].parent()
    V = RR^(r+1)
    UL = [V(UL[i]) for i in range(r)]
    return Matrix(RR, r, [u.dot_product(v) for u in UL for v in UL])

def _perp(alpha):
    """
    Return the orthogonal projection of alpha onto the trace zero space (i.e.
    onto the orthogonal complement of 1 under the Minkowski inner product)
    scaled by the degree of the number field in which alpha lives.
    """
    return alpha.parent().absolute_degree() * alpha - alpha.trace()

def BQF_to_H(Q):
    """
    Return the point in the upper-half plane corresponding to the positive-definite
    rational binary quadratic from Q. The point is returned as an element of QQbar.
    """
    G = Q.Gram_matrix_rational()
    xx = G[0,1] / G[0,0]
    yy = (G[1,1] / G[0,0] - xx^2).sqrt()
    return xx + QQbar(-1).sqrt() * yy
