# ---------------------------------------------------------------------------
# Phase 1.3: Transfer-Matrix BGFs, Graph Walks, Dynamical Zeta & Random Mappings
# ---------------------------------------------------------------------------

function _is_row_stochastic(T::AbstractMatrix{<:Real}; atol::Real=1e-6)
    m, n = size(T)
    m == n || return false
    for i in 1:m
        s = 0.0
        for j in 1:m
            T[i, j] >= -atol || return false
            s += T[i, j]
        end
        abs(s - 1.0) <= atol || return false
    end
    return true
end

function _stationary_dist(T::AbstractMatrix{<:Real})
    m = size(T, 1)
    M = Matrix{Float64}(I - T')
    M[end, :] .= 1.0
    rhs = zeros(Float64, m)
    rhs[end] = 1.0
    if abs(det(M)) > 1e-11
        pi_vec = M \ rhs
        if all(x -> x >= -1e-9, pi_vec)
            pi_vec = max.(pi_vec, 0.0)
            s = sum(pi_vec)
            s > 0 && return pi_vec ./ s
        end
    end
    # Fallback: Cesaro average of power iteration (handles reducible & periodic chains)
    v = fill(1.0 / m, m)
    acc = copy(v)
    Tt = Matrix{Float64}(T')
    steps = 250
    for _ in 1:steps
        v = Tt * v
        acc .+= v
    end
    acc ./= sum(acc)
    return acc
end

function _matmul_generic(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    r, k = size(A)
    _, c = size(B)
    C = Matrix{T}(undef, r, c)
    @inbounds for j in 1:c, i in 1:r
        s = zero(T)
        for p in 1:k
            s += A[i, p] * B[p, j]
        end
        C[i, j] = s
    end
    return C
end

function _matvec_generic(A::AbstractMatrix{T}, v::AbstractVector{T}) where {T}
    r, c = size(A)
    y = Vector{T}(undef, r)
    @inbounds for i in 1:r
        s = zero(T)
        for j in 1:c
            s += A[i, j] * v[j]
        end
        y[i] = s
    end
    return y
end

function _matpow_generic(A::AbstractMatrix{T}, n::Integer) where {T}
    m = size(A, 1)
    res = Matrix{T}(undef, m, m)
    for j in 1:m, i in 1:m
        res[i, j] = i == j ? one(A[1, 1]) : zero(A[1, 1])
    end
    base = copy(A)
    k = Int(n)
    while k > 0
        if isodd(k)
            res = _matmul_generic(res, base)
        end
        k >>= 1
        if k > 0
            base = _matmul_generic(base, base)
        end
    end
    return res
end


"""
    transfer_matrix_moments(T::AbstractMatrix{<:Real}, mark_matrix::AbstractMatrix{<:Real}, n::Integer;
                            init_dist=nothing, final_vec=nothing)

Compute exact analytical moments `(mean, var, std, count)` for the cumulative marked transition
statistic ``\\chi = \\sum_{t=1}^n M_{s_{t-1}, s_t}`` over paths of `n` transitions on an
``m \\times m`` transfer/transition matrix `T` (Flajolet & Sedgewick, Theorem V.7 & Section V.6):

``G(z, u) = \\mathbf{p}_0^T \\bigl(\\mathbf{I} - z \\mathbf{T}(u)\\bigr)^{-1} \\mathbf{f}, \\qquad T_{ij}(u) = T_{ij} \\, u^{M_{ij}}.``

Uses exact degree-2 `Taylor1` propagation in ``u = 1 + du`` in ``O(\\min(n m^2, m^3 \\log n))`` time.
If `init_dist` is `nothing`, defaults to the stationary distribution when `T` is row-stochastic,
and `ones(m)` otherwise. `final_vec` defaults to `ones(m)`.
"""
function transfer_matrix_moments(
    T::AbstractMatrix{<:Real},
    mark_matrix::AbstractMatrix{<:Real},
    n::Integer;
    init_dist=nothing,
    final_vec=nothing
)
    m, m2 = size(T)
    m == m2 || throw(ArgumentError("Transfer matrix T must be square, got $(size(T))"))
    size(mark_matrix) == (m, m) || throw(ArgumentError("mark_matrix size $(size(mark_matrix)) must match T ($m, $m)"))
    n >= 0 || throw(ArgumentError("Path length n must be >= 0, got $n"))

    p0 = if init_dist === nothing
        _is_row_stochastic(T) ? _stationary_dist(T) : ones(Float64, m)
    else
        length(init_dist) == m || throw(ArgumentError("init_dist length must be $m"))
        Float64.(init_dist)
    end
    f = if final_vec === nothing
        ones(Float64, m)
    else
        length(final_vec) == m || throw(ArgumentError("final_vec length must be $m"))
        Float64.(final_vec)
    end

    if n == 0
        cnt = Float64(dot(p0, f))
        return (mean = 0.0, var = 0.0, std = 0.0, count = cnt)
    end

    du = TaylorSeries.Taylor1(Float64, 2)
    u = 1.0 + du
    Tu = Matrix{TaylorSeries.Taylor1{Float64}}(undef, m, m)
    for j in 1:m, i in 1:m
        tij = Float64(T[i, j])
        mij = Float64(mark_matrix[i, j])
        if tij == 0.0
            Tu[i, j] = zero(du)
        elseif mij == 0.0
            Tu[i, j] = TaylorSeries.Taylor1(tij, 2)
        elseif mij == 1.0
            Tu[i, j] = tij * u
        else
            Tu[i, j] = tij * (u^mij)
        end
    end

    w = TaylorSeries.Taylor1{Float64}[TaylorSeries.Taylor1(fi, 2) for fi in f]
    if n <= 64 || n <= 2 * m
        for _ in 1:n
            w = _matvec_generic(Tu, w)
        end
    else
        P = _matpow_generic(Tu, n)
        w = _matvec_generic(P, w)
    end

    val = zero(du)
    for i in 1:m
        val += p0[i] * w[i]
    end

    a_n = Float64(TaylorSeries.getcoeff(val, 0))
    f1 = Float64(TaylorSeries.getcoeff(val, 1)) / a_n
    f2 = 2.0 * Float64(TaylorSeries.getcoeff(val, 2)) / a_n
    var_val = max(0.0, f2 + f1 - f1^2)
    return (mean = f1, var = var_val, std = sqrt(var_val), count = a_n)
end


function _solve_symbolic_linear(A_in::Matrix{SymPy.Sym}, b_in::Vector{SymPy.Sym})
    m = size(A_in, 1)
    A = copy(A_in)
    b = copy(b_in)
    for col in 1:m
        piv = 0
        for r in col:m
            if !isequal(SymPy.simplify(A[r, col]), SymPy.Sym(0))
                piv = r
                break
            end
        end
        piv == 0 && error("Singular matrix in transfer_matrix_gf")
        if piv != col
            A[col, :], A[piv, :] = A[piv, :], A[col, :]
            b[col], b[piv] = b[piv], b[col]
        end
        diag_el = A[col, col]
        for r in (col + 1):m
            if !isequal(A[r, col], SymPy.Sym(0))
                factor = A[r, col] / diag_el
                for c in col:m
                    A[r, c] = A[r, c] - factor * A[col, c]
                end
                b[r] = b[r] - factor * b[col]
            end
        end
    end
    x = Vector{SymPy.Sym}(undef, m)
    for i in m:-1:1
        s = b[i]
        for j in (i + 1):m
            s = s - A[i, j] * x[j]
        end
        x[i] = SymPy.simplify(s / A[i, i])
    end
    return x
end


"""
    transfer_matrix_gf(T::AbstractMatrix{<:Real}, z;
                       mark_matrix=nothing, init_dist=nothing, final_vec=nothing, u=1)

Evaluate the Transfer-Matrix Generating Function
``G(z, u) = \\mathbf{p}_0^T \\bigl(\\mathbf{I} - z \\mathbf{T}(u)\\bigr)^{-1} \\mathbf{f}``
for symbolic `SymPy.Sym` arguments, `TaylorSeries.Taylor1` series expansion, or numeric scalars.
"""
function transfer_matrix_gf(
    T::AbstractMatrix{<:Real},
    z;
    mark_matrix=nothing,
    init_dist=nothing,
    final_vec=nothing,
    u=1
)
    m, m2 = size(T)
    m == m2 || throw(ArgumentError("Transfer matrix T must be square, got $(size(T))"))

    p0 = if init_dist === nothing
        _is_row_stochastic(T) ? _stationary_dist(T) : ones(Float64, m)
    else
        Float64.(init_dist)
    end
    f = final_vec === nothing ? ones(Float64, m) : Float64.(final_vec)

    if z isa SymPy.Sym || u isa SymPy.Sym
        z_s = SymPy.Sym(z)
        u_s = SymPy.Sym(u)
        A = Matrix{SymPy.Sym}(undef, m, m)
        for i in 1:m, j in 1:m
            t_ij = SymPy.Sym(T[i, j])
            if mark_matrix !== nothing && mark_matrix[i, j] != 0
                t_ij = t_ij * (u_s ^ mark_matrix[i, j])
            end
            A[i, j] = (i == j ? SymPy.Sym(1) : SymPy.Sym(0)) - z_s * t_ij
        end
        sol = _solve_symbolic_linear(A, SymPy.Sym.(f))
        res = sum(SymPy.Sym(p0[i]) * sol[i] for i in 1:m)
        return SymPy.simplify(res)
    elseif z isa TaylorSeries.Taylor1
        ord = z.order
        Tu = Matrix{typeof(float(u) * 1.0)}(undef, m, m)
        for i in 1:m, j in 1:m
            tij = Float64(T[i, j])
            if mark_matrix !== nothing && mark_matrix[i, j] != 0
                Tu[i, j] = tij * (u ^ mark_matrix[i, j])
            else
                Tu[i, j] = tij * one(u)
            end
        end
        w = typeof(one(u) * 1.0).(f)
        acc = zero(z) + dot(p0, w)
        zk = z
        for k in 1:ord
            w = _matvec_generic(Tu, w)
            acc = acc + dot(p0, w) * zk
            if k < ord
                zk = zk * z
            end
        end
        return acc
    else
        Tu = Matrix{promote_type(Float64, typeof(u))}(undef, m, m)
        for i in 1:m, j in 1:m
            tij = Float64(T[i, j])
            if mark_matrix !== nothing && mark_matrix[i, j] != 0
                Tu[i, j] = tij * (u ^ mark_matrix[i, j])
            else
                Tu[i, j] = tij
            end
        end
        M = Matrix{eltype(Tu)}(I, m, m) - z * Tu
        return dot(p0, M \ f)
    end
end


"""
    graph_walks_count(A::AbstractMatrix{<:Real}, n::Integer; source=nothing, target=nothing)

Count walks of exact length `n` on a graph with adjacency or weight matrix `A`
(Flajolet & Sedgewick, Section V.5). Optionally restrict to walks starting at `source`
and/or ending at `target`.
"""
function graph_walks_count(A::AbstractMatrix{<:Real}, n::Integer; source=nothing, target=nothing)
    m, m2 = size(A)
    m == m2 || throw(ArgumentError("Adjacency matrix A must be square, got $(size(A))"))
    n >= 0 || throw(ArgumentError("Walk length n must be >= 0, got $n"))
    An = _matpow_generic(Matrix{promote_type(eltype(A), Int)}(A), n)
    if source !== nothing && target !== nothing
        return An[source, target]
    elseif source !== nothing
        return sum(@view An[source, :])
    elseif target !== nothing
        return sum(@view An[:, target])
    else
        return sum(An)
    end
end


"""
    graph_closed_walks_gf(A::AbstractMatrix{<:Real}, max_n::Integer)

Analyze closed walks, dynamical zeta polynomial, and topological entropy on a directed
or undirected graph with adjacency/transition matrix `A` (Flajolet & Sedgewick, Section V.5.1):

Returns a `NamedTuple`:
- `closed_walks`: vector `[Tr(A^1), Tr(A^2), ..., Tr(A^max_n)]` (coefficients of ``z^k`` in ``\\text{Tr}(zA(I - zA)^{-1})``)
- `zeta_inv_coeffs`: coefficients `[1, c_1, ..., c_m]` of the reciprocal dynamical zeta polynomial
  ``\\zeta_A(z)^{-1} = \\det(I - z A) = 1 + c_1 z + \\dots + c_m z^m`` via Newton's identities
- `spectral_radius`: dominant eigenvalue modulus ``\\rho(A)``
- `topological_entropy`: ``h_{\\text{top}} = \\log \\rho(A)``
"""
function graph_closed_walks_gf(A::AbstractMatrix{<:Real}, max_n::Integer)
    m, m2 = size(A)
    m == m2 || throw(ArgumentError("Matrix A must be square, got $(size(A))"))
    max_n >= 1 || throw(ArgumentError("max_n must be >= 1, got $max_n"))

    Af = Float64.(A)
    need_k = max(Int(max_n), m)
    traces = Vector{Float64}(undef, need_k)
    Ak = copy(Af)
    for k in 1:need_k
        traces[k] = tr(Ak)
        if k < need_k
            Ak = Ak * Af
        end
    end

    # Newton's identities for det(I - z A) = 1 + c_1 z + ... + c_m z^m
    c = zeros(Float64, m + 1)
    c[1] = 1.0
    for k in 1:m
        s = 0.0
        for j in 1:k
            s += c[k - j + 1] * traces[j]
        end
        c[k + 1] = -s / k
    end

    rho = maximum(abs.(eigvals(Af)))
    h_top = rho > 0 ? log(rho) : -Inf

    return (
        closed_walks = traces[1:Int(max_n)],
        zeta_inv_coeffs = c,
        spectral_radius = rho,
        topological_entropy = h_top
    )
end


"""
    random_mapping_moments(n::Integer)

Exact finite-``n`` and asymptotic statistics for a uniformly random discrete dynamical mapping
``f: \\{1, \\dots, n\\} \\to \\{1, \\dots, n\\}`` (Flajolet & Sedgewick, Section II.5.2 & Section VII.3.3):

Returns a `NamedTuple`:
- `expected_cyclic_nodes`: exact ``\\mathbb{E}_n[\\text{periodic states}] = Q(n) = \\sum_{k=1}^n \\frac{n!}{(n-k)! n^k} \\sim \\sqrt{\\pi n / 2}``
- `expected_components`: exact ``\\mathbb{E}_n[\\text{attractors}] = \\sum_{k=1}^n \\frac{1}{k} \\frac{n!}{(n-k)! n^k} \\sim \\frac{1}{2}(\\log(2n) + \\gamma)``
- `prob_connected`: exact probability of a single connected attractor ``Q(n) / n \\sim \\sqrt{\\pi / (2n)}``
- `expected_leaves`: exact expected number of transient terminal states (in-degree 0) ``n (1 - 1/n)^n \\sim n / e``
- `expected_fixed_points`: exact ``1.0``
- `cyclic_nodes_asym`: ``\\sqrt{\\pi n / 2}``
- `components_asym`: ``0.5 (\\log(2n) + \\gamma)``
- `tail_len_asym`: ``\\sqrt{\\pi n / 8}`` (mean transient distance ``\\lambda`` to limit cycle)
- `cycle_len_asym`: ``\\sqrt{\\pi n / 8}`` (mean limit cycle length ``\\mu`` from random state)
- `rho_len_asym`: ``\\sqrt{\\pi n / 2}`` (mean total orbit length ``\\rho = \\lambda + \\mu``)
- `largest_component_asym`: ``0.7582 n``
- `largest_tree_asym`: ``0.48 n``
"""
function random_mapping_moments(n::Integer)
    n >= 1 || throw(ArgumentError("State space size n must be >= 1, got $n"))
    nf = Float64(n)
    t_k = 1.0
    q_sum = 0.0
    comp_sum = 0.0
    for k in 1:Int(n)
        if k > 1
            t_k *= (nf - k + 1.0) / nf
        end
        q_sum += t_k
        comp_sum += t_k / k
        if t_k < 1e-16
            break
        end
    end

    euler_gamma = 0.5772156649015329
    return (
        expected_cyclic_nodes = q_sum,
        expected_components = comp_sum,
        prob_connected = q_sum / nf,
        expected_leaves = nf * ((1.0 - 1.0 / nf)^Int(n)),
        expected_fixed_points = 1.0,
        cyclic_nodes_asym = sqrt(π * nf / 2.0),
        components_asym = 0.5 * (log(2.0 * nf) + euler_gamma),
        tail_len_asym = sqrt(π * nf / 8.0),
        cycle_len_asym = sqrt(π * nf / 8.0),
        rho_len_asym = sqrt(π * nf / 2.0),
        largest_component_asym = 0.7582 * nf,
        largest_tree_asym = 0.48 * nf
    )
end
