def shape_of_a_number_field(K, check = True, prec = 50):
    """
    Return the shape of the number field.
    
    The shape of a number field of degree d is by definition the equivalence class of the (d-1)-dimensional lattice given by the orthogonal projection of the ring of integers of K onto the trace zero space (i.e. onto the orthogonal complement of 1) (equivalence means up to change of basis, orthogonal transformation, and homothety). The output of this function is the Gram matrix of the lattice. By default, this function checks if K is totally real or CM, in which cases it returns the Gram matrix over the rational numbers. Otherwise (or if check is set to False), a numerical approximation of the Gram matrix in RR is returned. In the latter case, the parameter prec is used for the precision of the real field. In order to provide a non-numerical answer in the general case, one would have to construct the normal closure (or at least the closure under complex conjugations). This is not practical at this time.
    
    """
    
    if K.absolute_degree() == 1:
        return matrix(QQ, 0)
    
    if check:
        if K.is_totally_real():
            return shape_of_a_totally_real_number_field(K)
        try:
            if K.is_CM():
                return shape_of_a_CM_number_field(K)
        except(AttributeError):    #temporary: only needed if trac #11770 has not been applied
            if CM_field_functionality.is_CM(K):
                return shape_of_a_CM_number_field(K)
    
    #Otherwise do it with AA and QQbar
    if K.is_relative():
        KK = K.absolute_field(K.variable_names()[0] * 2)
    else:
        KK = K
    OK = KK.maximal_order()
    B = OK.basis()
    if B[0] != 1:
        B = _sub_one_into_basis(OK)
    Bperp = []
    for a in range(1, KK.absolute_degree()):
        Bperp.append(_perp(B[a]))
    sigmas = list(K.embeddings(AA))
    iotas = list(K.embeddings(QQbar))
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
    
    OK = K.maximal_order()
    B = OK.basis()
    if B[0] != 1:
        B = _sub_one_into_basis(OK)
    Bperp = []
    for a in range(1, K.absolute_degree()):
        Bperp.append(_perp(B[a]))
    G = matrix(QQ, K.absolute_degree() - 1, K.absolute_degree() - 1)
    for a in range(len(Bperp)):
        for j in range(a, len(Bperp)):
            G[a,j] = (Bperp[a] * Bperp[j]).trace()
            if a != j:
                G[j,a] = G[a,j]
    return G

def minkowski_vector(alpha, sigmas, taus):
    n = len(sigmas) + 2 * len(taus)
    V = AA^n
    return V([sigma(alpha) for sigma in sigmas] + sum([[tau(alpha).real(), tau(alpha).imag()] for tau in taus], []))

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

def _sub_one_into_basis(OK):
    """
    The integral basis function does not always return a basis containing
    1. This function receives such a basis and modifies it by substituting
    1 for one of the basis elements. The returned basis contains 1 as its
    first element.
    
    """
    
    B = []
    for v in OK.basis():
        B.append(v)
    w = OK.coordinates(OK.one())
    sub_index = None
    for i in range(len(w)):
        if w[i] == 1 or w[i] == -1:
            sub_index = i
            break
    if not sub_index is None:
        B[sub_index] = B[0]
        B[0] = OK.one()
        #Debug test: print K.discriminant() == K.discriminant(B)
        return B
    
    #Otherwise, more complicated
    found = False
    for i in range(len(w) - 1):
        for j in range(i + 1, len(w)):
            g, x, y = xgcd(w[i], w[j])
            if g == 1:
                found = True
                break
        if found:
            break
    if not found:    #This should never happen
        raise ValueError("You've given me a field whose integral basis \
                cannot contain 1!")
    
    companion = -y * B[i] + x * B[j]
    B[i] = B[0]
    B[j] = companion
    B[0] = OK.one()
    #Debug test: print K.discriminant() == K.discriminant(B)
    return B
