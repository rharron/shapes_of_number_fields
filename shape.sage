from sage.rings.number_field.number_field_base import is_NumberField
from sage.libs.pari.convert_sage import gen_to_sage

def shape_of_a_number_field(K, check = True):
    """
    Return the shape of the number field.

    The shape of a number field of degree d is by definition the equivalence
    class of the (d-1)-dimensional lattice given by the orthogonal projection of
    the ring of integers of K onto the trace zero space (i.e. onto the orthogonal complement of 1) (equivalence means up to change of basis, orthogonal transformation, and homothety). The output of this function is the Gram matrix of the lattice. By default, this function checks if K is totally real in which case it returns the Gram matrix over the integers. Otherwise (or if check is set to False), the algebraic reals are used.

    """

    if K.absolute_degree() == 1:
        return matrix(ZZ, 0)

    if K.absolute_degree() == 2:
        return Matrix(ZZ, 1, [1])

    if check:
        if K.is_totally_real():
            return shape_of_a_totally_real_number_field(K)
        #try:
        #    if K.is_CM():
        #        return shape_of_a_CM_number_field(K)
        #except(AttributeError):    #temporary: only needed if trac #11770 has not been applied
        #    if CM_field_functionality.is_CM(K):
        #        return shape_of_a_CM_number_field(K)

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
    #ME = KK.Minkowski_embedding(B = [1] + Bperp, prec = prec)
    #return ME.submatrix(0,1,-1,-1).transpose() * ME.submatrix(0,1,-1,-1)

def shape_of_a_totally_real_number_field(K, check = False):
    """
    Return the shape of the totally real number field K.

    The shape of a totally real number field is a quadratic form over the rational numbers and the output of this function is the Gram matrix (over QQ). The computations are purely algebraic as the Minkowski inner product is simply the trace form. By default, this function does not check that the given field is totally real. To force a check, set "check" to True.

    """

    if check:
        if not K.is_totally_real():
            raise ValueError("Field passed to shape_of_a_totally_real_number_field is not totally real.")
    return trace_zero_form(K)

def trace_zero_form(K):
    r"""

    INPUT:

    - ``K`` -- a number field or an irreducible polynomial over the rationals

    OUTPUT:

    - The Gram matrix projection of the trace zero form to the trace zero space with the definition
    that always gives an integral quadratic form.
    """
    if is_NumberField(K):
        T2 = K.pari_nf()[4][3]
    else: #assumes K is a polynomial over QQ or I guess anything that can be converted to pari and passed to nfinit
        K = pari(K)
        T2 = K.nfinit()[4][3]
    G = gen_to_sage(T2)
    return _gram_to_perp_gram(G)

def unit_shape(K, algorithm='pari'):
    """
    pari is much, much faster. the other option for algorithm is 'sage'
    """
    if algorithm == 'pari':
        return unit_shape_pari(K)
    elif algorith == 'sage':
        return unit_shape_sage(K)
    raise ValueError("Parameter algorithm in unit_shape must be either \'pari\' or \'sage\'.")

def _gram_to_perp_gram(G):
    """
    Input an nxn Gram matrix where the first row and column are assumed to correspond to the number 1
    and output an (n-1)x(n-1) matrix that is the Gram matrix of the perps of the remaining elements.
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
    n = len(sigmas) + 2 * len(taus)
    V = AA^n
    return V([sigma(alpha) for sigma in sigmas] + sum([[tau(alpha).real(), tau(alpha).imag()] for tau in taus], []))

def multiplicative_minkowski_vector(u, sigmas, taus):
    m = len(sigmas) + len(taus)
    v = [sigma(u).abs().log() for sigma in sigmas] + [(tau(u).abs()^2).log() for tau in taus]
    return (RR^m)(v)

def unit_shape_sage(K):
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
    Return the orthogonal projection of alpha onto the trace zero space (i.e. onto the orthogonal complement of 1 under the Minkowski inner product).
    """
    return alpha.parent().absolute_degree() * alpha - alpha.trace()

def BQF_to_H(Q):
    G = Q.Gram_matrix_rational()
    xx = G[0,1] / G[0,0]
    yy = (G[1,1] / G[0,0] - xx^2).sqrt()
    return xx + QQbar(-1).sqrt() * yy
