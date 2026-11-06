// Optional categorical-FE Schur projection. Eliminate the largest block,
// then solve B' W M_A B using diagonally preconditioned conjugate gradients.
// No dense dummy or Schur matrix is formed. The two-effect specialization
// applies its operator to coalesced cells; the multiway path scans observations.
mata:
class PpmlSchurGroup
{
    real colvector levels, order, cell_order, cell_a, cell_b, ca_order, cb_order
    real matrix info, cell_info, ca_info, cb_info
    real scalar first, last
}
class PpmlSchur
{
    real scalar enabled, eliminated, n, size, iterations
    real scalar calls, fallbacks, cg_iterations
    real colvector a_levels, a_order, a_weight, weight, roots, diagonal, cells
    real matrix a_info
    class PpmlSchurGroup rowvector groups
    void init(), project()
    real matrix within(), transpose(), multiply()
    real scalar partial_out()
}

void PpmlSchur::init(class FixedEffects scalar H)
{
    class Factor scalar A, B, C, CA, CB
    real scalar g, j
    enabled = iterations = calls = fallbacks = cg_iterations = 0
    if (H.G<2 | any(H.num_slopes) | !all(H.has_intercept) | H.prune) return
    // RRE compares a particular original-input iteration history.
    if (H.compute_rre | H.storing_alphas) return
    eliminated = 1
    for (g=2; g<=H.G; g++) {
        if (H.factors[g].num_levels>H.factors[eliminated].num_levels) eliminated=g
    }
    A = H.factors[eliminated]
    n = A.num_obs
    a_levels = A.levels
    a_order = A.sort((1::n))
    a_info = A.info
    groups = PpmlSchurGroup(H.G-1)
    size = 0
    j = 0
    for (g=1; g<=H.G; g++) {
        if (g==eliminated) continue
        ++j
        B = H.factors[g]
        groups[j].levels = B.levels
        groups[j].order = B.sort((1::n))
        groups[j].info = B.info
        groups[j].first = size+1
        size = size+B.num_levels
        groups[j].last = size
        C = join_factors(A, B, 1, 1, 1)
        C.panelsetup()
        groups[j].cell_order = C.sort((1::n))
        groups[j].cell_info = C.info
        groups[j].cell_a = C.keys[.,1]
        groups[j].cell_b = C.keys[.,2]
        CB = _factor(C.keys[.,2], 1, 0, "", 1, 1, ., 0)
        CB.panelsetup()
        groups[j].cb_order = CB.sort((1::C.num_levels))
        groups[j].cb_info = CB.info
        // Only the two-effect operator needs to sum cells within A.
        if (H.G==2) {
            CA = _factor(C.keys[.,1], 1, 0, "", 1, 1, ., 0)
            CA.panelsetup()
            groups[j].ca_order = CA.sort((1::C.num_levels))
            groups[j].ca_info = CA.info
        }
    }
    enabled = 1
}

real matrix PpmlSchur::within(real matrix x)
{
    return(x - (panelsum((weight:*x)[a_order,.],a_info):/a_weight)[a_levels,.])
}

real matrix PpmlSchur::transpose(real matrix x)
{
    real matrix out
    real scalar g
    out = J(size,cols(x),0)
    for (g=1; g<=cols(groups); g++) {
        out[|groups[g].first,1 \ groups[g].last,cols(x)|] =
            panelsum((weight:*x)[groups[g].order,.],groups[g].info)
    }
    return(out:/roots)
}

real matrix PpmlSchur::multiply(real matrix v)
{
    real matrix x, coefficients
    real scalar g
    if (cols(groups)==1) {
        x = panelsum((cells:*v[groups[1].cell_b,.])[groups[1].ca_order,.],groups[1].ca_info)
        return(v-panelsum((cells:*x[groups[1].cell_a,.])[groups[1].cb_order,.],groups[1].cb_info))
    }
    coefficients = v:/roots
    x = J(n,cols(v),0)
    for (g=1; g<=cols(groups); g++) {
        x = x+coefficients[groups[g].first:-1:+groups[g].levels,.]
    }
    return(transpose(within(x)))
}

// Transactional kernel: failures leave data untouched. Weights and columns
// are rescaled internally, without changing their projection, to limit overflow.
real scalar PpmlSchur::partial_out(real matrix data, real colvector w,
                                 real scalar tolerance, real scalar maxiter)
{
    real scalar g, i
    real colvector cw, nb, fraction, invdiag, idx
    real rowvector scale, norm0, target, rho, previous_rho, alpha, beta, active, denom
    real matrix x, b, a, r, z, p, ap, candidate, coefficients, moment
    iterations = 0
    if (!enabled | rows(data)!=n | !cols(data) | hasmissing(data) | hasmissing(w)) return(0)
    if (rows(w)!=1 & rows(w)!=n) return(0)
    if (any(w:<=0) | missing(tolerance) | tolerance<=0 | missing(maxiter) | maxiter<0) return(0)
    weight = rows(w)==1 ? J(n,1,1) : w:/max(w)
    if (any(weight:<=0)) return(0) // underflow across an extreme weight range
    a_weight = panelsum(weight[a_order],a_info)
    roots = diagonal = J(size,1,0)
    for (g=1; g<=cols(groups); g++) {
        idx = (groups[g].first::groups[g].last)
        nb = panelsum(weight[groups[g].order],groups[g].info)
        roots[idx] = sqrt(nb)
        cw = panelsum(weight[groups[g].cell_order],groups[g].cell_info)
        fraction = cw:/a_weight[groups[g].cell_a]
        if (any(fraction:>1+64*epsilon(1))) return(0)
        fraction = rowmin((fraction,J(rows(fraction),1,1)))
        diagonal[idx] = panelsum((cw:*(1:-fraction))[groups[g].cb_order],groups[g].cb_info):/nb
        if (cols(groups)==1) cells = (cw:/sqrt(a_weight[groups[g].cell_a])):/sqrt(nb[groups[g].cell_b])
    }
    if (hasmissing(roots) | any(roots:<=0) | hasmissing(diagonal) | any(diagonal:<0)) return(0)
    invdiag = J(size,1,0)
    idx = selectindex(diagonal:>0)
    if (rows(idx)) invdiag[idx] = 1:/diagonal[idx]
    if (hasmissing(invdiag)) return(0)
    scale = colmax(abs(data))
    scale = scale+(scale:==0)
    x = data:/scale
    norm0 = sqrt(colsum(weight:*x:^2))
    target = tolerance:*rowmax((norm0',J(cols(data),1,1e-30)))'
    b = transpose(within(x))
    a = J(size,cols(data),0)
    r = b
    // A zero diagonal represents a redundant level. Its RHS must be zero too.
    idx = selectindex(diagonal:==0)
    if (rows(idx)) {
        if (any(sqrt(colsum(r[idx,.]:^2)):>target)) return(0)
        r[idx,.] = J(rows(idx),cols(data),0)
    }
    active = sqrt(colsum(r:^2)):>target
    z = invdiag:*r
    p = z:*active
    rho = colsum(r:*z)
    for (i=1; i<=maxiter & any(active); i++) {
        ap = multiply(p)
        denom = colsum(p:*ap)
        if (hasmissing(denom) | any((denom:<=0):&active)) return(0)
        alpha = editmissing(rho:/denom,0):*active
        a = a+p:*alpha
        r = r-ap:*alpha
        iterations = i
        active = sqrt(colsum(r:^2)):>target
        if (!any(active)) break
        z = invdiag:*r
        previous_rho = rho
        rho = colsum(r:*z)
        beta = editmissing(rho:/previous_rho,0)
        p = (z+p:*beta):*active
    }
    if (any(active) | hasmissing(a)) return(0)
    // Check a freshly evaluated reduced residual, not just the CG recurrence.
    if (any(sqrt(colsum((b-multiply(a)):^2)):>10:*target)) return(0)
    coefficients = a:/roots
    candidate = x
    for (g=1; g<=cols(groups); g++) {
        candidate = candidate-coefficients[groups[g].first:-1:+groups[g].levels,.]
    }
    candidate = within(candidate)
    if (hasmissing(candidate)) return(0)
    moment = transpose(candidate)
    if (any(sqrt(colsum(moment:^2)):>10:*target)) return(0)
    moment = panelsum((weight:*candidate)[a_order,.],a_info):/sqrt(a_weight)
    if (any(sqrt(colsum(moment:^2)):>10:*target)) return(0)
    if (any(sqrt(colsum(weight:*candidate:^2)):>norm0:+10:*target)) return(0)
    candidate = candidate:*scale
    if (hasmissing(candidate)) return(0)
    swap(data,candidate)
    return(1)
}

// Retain reghdfe's original-input norm, means, and zero-column convention.
// Both call sites use save_tss=standardize=first_is_depvar=0, flush=1.
// Near either FE- or regressor-collinearity thresholds, use the original path.
void PpmlSchur::project(class FixedEffects scalar H, real matrix data, real colvector w)
{
    real matrix candidate, gram
    real rowvector norm0, norm1, ratio, keep, zero, means, scales
    real scalar cutoff, ok
    if (!enabled) {
        H._partial_out(data,0,0,0,1)
        return
    }
    ++calls
    candidate = data
    norm0 = diagonal(cross(data,data))'
    means = H.compute_constant ? mean(data,H.weight) : J(1,cols(data),1)
    cutoff = min((1e-6,H.tolerance/10))
    ok = partial_out(candidate,w,H.tolerance,H.maxiter)
    cg_iterations = cg_iterations+iterations
    if (ok) {
        norm1 = diagonal(cross(candidate,candidate))'
        ratio = norm1:/norm0
        ok = !hasmissing(norm0) & !hasmissing(norm1) & all(norm0:>0)
        ok = ok & !any((ratio:>=cutoff/100):&(ratio:<=100*cutoff))
    }
    if (ok) {
        H._partial_out(candidate,0,0,0,1)
        ratio = diagonal(cross(candidate,candidate))':/norm0
        ok = !any((ratio:>=cutoff/100):&(ratio:<=100*cutoff))
        keep = ratio:>cutoff
        zero = selectindex(!keep)
        if (cols(zero)) candidate[.,zero] = J(rows(data),cols(zero),0)
        // A small K-by-K condition check protects the separate invsym rank
        // decision. Keep the zero-column fallback unconditional: polishing
        // used candidate norms and could otherwise falsely zero a column
        // after a legitimate weighted increase in its unweighted norm.
        gram = quadcross(candidate,weight,candidate)
        scales = sqrt(diagonal(gram))'
        if (any(scales:<=0) | hasmissing(scales)) ok=0
        else {
            gram = (gram:/scales):/scales'
            ok = ok & cond(gram)<1e10
        }
    }
    if (!ok) {
        ++fallbacks
        H._partial_out(data,0,0,0,1)
        return
    }
    H.kept = keep
    H.means = means
    swap(data,candidate)
}
end
