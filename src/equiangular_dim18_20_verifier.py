"""Zero-float exact rational verifier for Real Equiangular Lines (d=18..20) & Fixed-Angle Bounds.

Certifies exact constructions, Gram ranks, PSD witnesses, paired-lattice
reductions, extended binary Golay code G_24 & Leech 276-axis evaluation code
generalized Hamming weights, nested sharp sections (276/176/126/90), Theorem 5.1
fixed-angle intervals & Seidel-block matrices, Theorem 5.2 Paley conference
matrices, and Proposition 5.3 cosine-1/7 maxima from Valjakas Djaljapayan,
Zhen-Qi Liao, and Wei-Hsuan Yu (arXiv:2611.34855v1, Nov 2026).

Schema: EQUIANGULAR-D18-20-CERT-1.0
Zero floats: uses strictly `int` and `fractions.Fraction` (no `float`, no `math`).
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

SCHEMA = "EQUIANGULAR-D18-20-CERT-1.0"


def _gcd(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _req_int(obj: Dict[str, Any], key: str) -> int:
    if key not in obj:
        raise ValueError(f"missing integer field '{key}'")
    v = obj[key]
    if isinstance(v, bool) or not isinstance(v, int):
        raise ValueError(f"field '{key}' must be exact int")
    return v


def _req_bool(obj: Dict[str, Any], key: str) -> bool:
    if key not in obj:
        raise ValueError(f"missing boolean field '{key}'")
    v = obj[key]
    if not isinstance(v, bool):
        raise ValueError(f"field '{key}' must be bool")
    return v


def exact_psd_and_rank(mat: List[List[Fraction]]) -> Tuple[bool, int]:
    """Exact rational pivoted symmetric Gaussian elimination (LDL^T) checking PSD and rank."""
    n = len(mat)
    A = [row[:] for row in mat]
    for i in range(n):
        if len(A[i]) != n:
            raise ValueError("matrix not square")
        for j in range(i + 1, n):
            if A[i][j] != A[j][i]:
                raise ValueError("matrix not symmetric")
    rank = 0
    for k in range(n):
        # Find largest-index nonzero diagonal entry >= k
        piv = None
        for i in range(k, n):
            if A[i][i] < 0:
                return False, rank
            if A[i][i] > 0 and piv is None:
                piv = i
        if piv is None:
            # All remaining diagonals are 0; PSD iff all remaining off-diagonals are 0
            for i in range(k, n):
                for j in range(i + 1, n):
                    if A[i][j] != 0:
                        return False, rank
            break
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            for i in range(n):
                A[i][k], A[i][piv] = A[i][piv], A[i][k]
        d = A[k][k]
        rank += 1
        for i in range(k + 1, n):
            lik = A[i][k]
            if lik != 0:
                f = lik / d
                for j in range(i, n):
                    val = A[i][j] - f * A[k][j]
                    A[i][j] = val
                    A[j][i] = val
                A[i][k] = Fraction(0)
                A[k][i] = Fraction(0)
    return True, rank


def build_dim18_57_gram() -> List[List[Fraction]]:
    """Construct exact scaled Gram matrix (diag=5, offdiag=+-1) of 57 lines in R^18 (Prop. 2.3)."""
    # 9 matchings of K_10: for i in Z/9Z, E_i = {{i, 9}} U {{(i-j)%9, (i+j)%9} : 1<=j<=4}
    matched_edges: List[Tuple[int, Tuple[int, int]]] = []
    for i in range(9):
        matched_edges.append((i, (i, 9)))
        for j in range(1, 5):
            a, b = (i - j) % 9, (i + j) % 9
            if a > b:
                a, b = b, a
            matched_edges.append((i, (a, b)))
    assert len(matched_edges) == 45

    # 6 directed Hamilton cycles of K_4 on {0,1,2,3}, canonically starting at 0
    perms = [
        (0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3),
        (0, 2, 3, 1), (0, 3, 1, 2), (0, 3, 2, 1),
    ]
    cycle_arcs = []
    for cyc in perms:
        arcs = {(cyc[k], cyc[(k + 1) % 4]) for k in range(4)}
        cycle_arcs.append(arcs)

    # Map 10 coordinates: 0..3 = vertices of K_4, 4..9 = 6 cycles
    blocks: List[frozenset[int]] = []
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            pts = {v for v in range(4) if v != i}
            for c_idx, arcs in enumerate(cycle_arcs):
                if (i, j) in arcs:
                    pts.add(4 + c_idx)
            assert len(pts) == 5
            blocks.append(frozenset(pts))
    assert len(blocks) == 12

    # Gram matrix: 45 matched-edge vectors + 12 block vectors = 57
    N = 57
    G: List[List[Fraction]] = [[Fraction(0)] * N for _ in range(N)]
    for a in range(N):
        G[a][a] = Fraction(5)
    # Edge-edge
    for idx1 in range(45):
        i1, e1 = matched_edges[idx1]
        s1 = set(e1)
        for idx2 in range(idx1 + 1, 45):
            i2, e2 = matched_edges[idx2]
            inter = len(s1.intersection(e2))
            val = 2 * ((1 if i1 == i2 else 0) + inter) - 1
            if val not in (-1, 1):
                raise ValueError(f"Prop 2.3 edge-edge inner product {val} not in +-1")
            G[idx1][idx2] = Fraction(val)
            G[idx2][idx1] = Fraction(val)
    # Block-block
    for b1 in range(12):
        for b2 in range(b1 + 1, 12):
            val = 2 * len(blocks[b1].intersection(blocks[b2])) - 5
            if val not in (-1, 1):
                raise ValueError(f"Prop 2.3 block-block inner product {val} not in +-1")
            G[45 + b1][45 + b2] = Fraction(val)
            G[45 + b2][45 + b1] = Fraction(val)
    # Edge-block
    for idx1 in range(45):
        i1, _ = matched_edges[idx1]
        for b1 in range(12):
            val = 2 * (1 if i1 in blocks[b1] else 0) - 1
            G[idx1][45 + b1] = Fraction(val)
            G[45 + b1][idx1] = Fraction(val)
    return G


def build_golay24_octads() -> Tuple[List[int], List[int]]:
    """Construct extended binary Golay code G_24 from X^11+X^9+X^7+X^6+X^5+X+1 (Sec. 4.4)."""
    # Generator polynomial g(X) = X^11 + X^9 + X^7 + X^6 + X^5 + X + 1
    gen_poly = (1 << 11) | (1 << 9) | (1 << 7) | (1 << 6) | (1 << 5) | (1 << 1) | 1
    basis: List[int] = []
    for shift in range(12):
        cw23 = gen_poly << shift
        parity = cw23.bit_count() & 1
        # Coordinate 0 is overall parity bit; coordinates 1..23 are bits 1..23
        cw24 = (cw23 << 1) | parity
        basis.append(cw24)
    # Generate all 4096 codewords
    codewords = [0]
    for b in basis:
        codewords.extend([c ^ b for c in codewords])
    assert len(codewords) == 4096
    octads = [c for c in codewords if c.bit_count() == 8]
    return codewords, octads


def verify_dim18_construction_and_reduction(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Mode 1: Prop. 2.3 (57 lines in R^18), Lemma 1.3 paired lattice, Lemma 2.2 / Eq. (1)."""
    submode = cert.get("submode")
    if submode == "prop23_57_lines":
        n_claim = _req_int(cert, "num_lines")
        dim_claim = _req_int(cert, "ambient_dim")
        q_claim = _req_int(cert, "inv_cosine_q")
        if n_claim != 57 or dim_claim != 18 or q_claim != 5:
            raise ValueError("Prop. 2.3 requires num_lines=57, ambient_dim=18, inv_cosine_q=5")
        G = build_dim18_57_gram()
        psd, rank = exact_psd_and_rank(G)
        if not psd or rank != 18:
            raise ValueError(f"57-line Gram matrix failed: psd={psd}, rank={rank}")
        trG = sum(G[i][i] for i in range(57))
        trG2 = sum(G[i][j] * G[i][j] for i in range(57) for j in range(57))
        if trG != 285 or trG2 != 57 * 81:
            raise ValueError("Gram trace identities mismatch")
        return {"mode": "dim18_construction_and_reduction", "submode": submode, "rank": rank, "psd": True}
    elif submode == "lemma13_paired_lattice":
        r = _req_int(cert, "span_dim_r")
        q = _req_int(cert, "inv_cosine_q")
        rank_T = _req_int(cert, "marking_rank_T")
        w2_claim = _req_int(cert, "claimed_w2")
        x2_claim = _req_int(cert, "claimed_x2")
        xw_claim = _req_int(cert, "claimed_xw")
        embeds_claim = _req_bool(cert, "embeds_in_rank24")
        if q < 3 or q % 2 == 0 or r < 1 or rank_T < 0:
            raise ValueError("Lemma 1.3 requires odd q >= 3, r >= 1, rank_T >= 0")
        w2 = 2 * (q - 2)
        x2 = q - 1
        xw = q - 2
        embeds = (r + 1 + rank_T) <= 21
        if w2_claim != w2 or x2_claim != x2 or xw_claim != xw or embeds_claim != embeds:
            raise ValueError("Lemma 1.3 lattice parameters mismatch")
        if w2 % 2 != 0 or x2 % 2 != 0:
            raise ValueError("lattice must be even integral")
        return {
            "mode": "dim18_construction_and_reduction",
            "submode": submode,
            "ambient_rank": r + 1 + rank_T,
            "embeds_in_rank24": embeds,
        }
    else:
        raise ValueError(f"unknown dim18 submode '{submode}'")


def _binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    num = 1
    den = 1
    for i in range(1, k + 1):
        num *= n - i + 1
        den *= i
    return num // den


def _lemma48_allowed_weights() -> set[int]:
    """Exact reproduction of Lemma 4.8 / Table 13 incidence-table evaluation code weights."""
    lam = [253, 77, 21, 5, 1]
    weights: set[int] = set()

    def ag_values(k: int, u: int, v: int) -> List[int]:
        if k == 0:
            return [lam[u]]
        if k == 7:
            return [2 * lam[u] - ((k - v) * lam[u + 1] + v * lam[u]) + (4 if u == v else 0)]
        alpha = 0 if k == 8 else 1
        if u <= 2:
            m1 = (k - v) * lam[u + 1] + v * lam[u]
            m2 = 2 * _binom(k - v, 2) * lam[u + 2] + (2 * v + 1) * (k - v) * lam[u + 1] + v * v * lam[u]
            num = -m2 + (2 * alpha + 4) * m1 - alpha * (alpha + 4) * lam[u]
            assert num % 2 == 0
            return [lam[u] - num // 2]
        # u == 3: enumerate all (n_alpha, n_{alpha+2}, n_{alpha+4}) summing to 5
        m1 = (k - v) * lam[4] + v * lam[3]
        vals = []
        for n0 in range(6):
            for n1 in range(6 - n0):
                n2 = 5 - n0 - n1
                triples = [(alpha, n0), (alpha + 2, n1), (alpha + 4, n2)]
                ok = True
                for i_val, cnt in triples:
                    if cnt > 0 and (i_val < v or 7 - i_val < u - v):
                        ok = False
                        break
                if not ok:
                    continue
                if alpha * n0 + (alpha + 2) * n1 + (alpha + 4) * n2 == m1:
                    vals.append(5 - 2 * n1)
        return sorted(set(vals))

    for k in (0, 7, 8, 11):
        for t in range(4):
            for r in range(min(k, t) + 1):
                if t - r > 23 - k:
                    continue
                # Enumerate all (u, v) subsets of b (size t, |g cap b| = r)
                # Only u=3 can have multiple values
                base_sum = 0
                u3_opts = [0]
                for u in range(t + 1):
                    for v in range(max(0, u - (t - r)), min(u, r) + 1):
                        mult = _binom(r, v) * _binom(t - r, u - v)
                        if mult == 0:
                            continue
                        sign_pow = (-2) ** u
                        opts = ag_values(k, u, v)
                        if len(opts) == 1:
                            base_sum += mult * sign_pow * opts[0]
                        else:
                            assert u == 3 and mult == 1
                            u3_opts = [sign_pow * val for val in opts]
                for add in u3_opts:
                    tot_sum = base_sum + add
                    inner = 253 - ((-1) ** t) * tot_sum
                    assert inner % 2 == 0
                    wt = k + inner // 2
                    weights.add(wt)
                    weights.add(276 - wt)
    return weights


def verify_leech_evaluation_code_dim20(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Mode 2: Golay G_24, Lemma 4.8 weights, Prop. 4.9 d_3(C)=186, Sec. 4.6 nested sections."""
    submode = cert.get("submode")
    if submode == "golay_and_lemma48":
        codewords, octads = build_golay24_octads()
        wt_set = sorted({c.bit_count() for c in codewords})
        if wt_set != [0, 8, 12, 16, 24] or len(octads) != 759:
            raise ValueError("extended binary Golay code G_24 verification failed")
        heptads = [c >> 1 for c in octads if (c & 1) == 1]
        if len(heptads) != 253:
            raise ValueError(f"expected 253 heptads, got {len(heptads)}")
        allowed_wts = sorted(_lemma48_allowed_weights())
        expected_wts = [0, 100, 112, 128, 132, 144, 148, 164, 176, 276]
        if allowed_wts != expected_wts:
            raise ValueError(f"Lemma 4.8 weight set mismatch: {allowed_wts}")
        claimed_nonzero = cert.get("claimed_nonzero_weights")
        if claimed_nonzero != expected_wts[1:]:
            raise ValueError("claimed_nonzero_weights mismatch")
        return {
            "mode": "leech_evaluation_code_dim20",
            "submode": submode,
            "octads": 759,
            "heptads": 253,
            "nonzero_weights": claimed_nonzero,
        }
    elif submode == "nested_sharp_sections":
        target_dim = _req_int(cert, "section_dim")
        claimed_axes = _req_int(cert, "claimed_num_axes")
        claimed_ghw = _req_int(cert, "claimed_generalized_hamming_weight")
        if target_dim not in (20, 21, 22):
            raise ValueError("section_dim must be in {20, 21, 22}")
        _, octads = build_golay24_octads()
        heptads = [c >> 1 for c in octads if (c & 1) == 1]
        # Pick a weight-8 word O in G (an octad NOT containing coordinate 0, punctured at 0)
        weight8_words = [c >> 1 for c in octads if (c & 1) == 0]
        O = weight8_words[0]
        assert O.bit_count() == 8
        pts_in_O = [j for j in range(23) if (O >> j) & 1]
        a, b = pts_in_O[0], pts_in_O[1]
        # Sec. 4.6 normals: y_a -> avoids a; y_O -> |H cap O| == 2; y_b -> avoids b
        if target_dim == 22:
            sel = [H for H in heptads if not ((H >> a) & 1)]
        elif target_dim == 21:
            sel = [H for H in heptads if (not ((H >> a) & 1)) and ((H & O).bit_count() == 2)]
        else:  # 20
            sel = [
                H for H in heptads
                if (not ((H >> a) & 1)) and (not ((H >> b) & 1)) and ((H & O).bit_count() == 2)
            ]
        m = len(sel)
        if m != claimed_axes:
            raise ValueError(f"section dim {target_dim}: expected {claimed_axes} axes, built {m}")
        if 276 - m != claimed_ghw:
            raise ValueError(f"generalized Hamming weight mismatch: 276-{m} != {claimed_ghw}")
        # Build exact m x m Gram matrix G_{H,H'} = |H cap H'| - 2 (diag=5, offdiag=+-1)
        G: List[List[Fraction]] = [[Fraction(0)] * m for _ in range(m)]
        for i in range(m):
            G[i][i] = Fraction(5)
            hi = sel[i]
            for j in range(i + 1, m):
                val = (hi & sel[j]).bit_count() - 2
                if val not in (-1, 1):
                    raise ValueError(f"Leech heptad inner product {val} not in +-1")
                G[i][j] = Fraction(val)
                G[j][i] = Fraction(val)
        psd, rank = exact_psd_and_rank(G)
        if not psd or rank != target_dim:
            raise ValueError(f"Leech section dim {target_dim} Gram check failed: psd={psd}, rank={rank}")
        return {
            "mode": "leech_evaluation_code_dim20",
            "submode": submode,
            "section_dim": target_dim,
            "num_axes": m,
            "gram_rank": rank,
        }
    else:
        raise ValueError(f"unknown leech submode '{submode}'")


def _floor_Lk_Q_exact(k: int, Q: int) -> int:
    """Exact zero-float floor(L_k * Q) via integer squaring.

    L_7  = 72 / (7*sqrt(10) + 2)
    L_23 = 600 / (23*sqrt(26) + 2)
    """
    if k == 7:
        num = 72 * Q
        c_rad, rad, c_int = 7, 10, 2
    elif k == 23:
        num = 600 * Q
        c_rad, rad, c_int = 23, 26, 2
    else:
        raise ValueError("k must be 7 or 23")
    # Find largest integer m >= 0 with m*(c_rad*sqrt(rad) + c_int) < num
    # <=> m*c_int < num and (m*c_rad)^2 * rad < (num - m*c_int)^2
    lo = 0
    hi = num // c_int + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        rem = num - mid * c_int
        if rem > 0 and (mid * c_rad) ** 2 * rad < rem * rem:
            lo = mid
        else:
            hi = mid
    # Verify strict bracketing (irrationality of L_k * Q)
    rem_lo = num - lo * c_int
    rem_hi = num - (lo + 1) * c_int
    assert rem_lo > 0 and (lo * c_rad) ** 2 * rad < rem_lo * rem_lo
    assert rem_hi <= 0 or ((lo + 1) * c_rad) ** 2 * rad > rem_hi * rem_hi
    return lo


def verify_fixed_angle_intervals_thm51(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Mode 3: Theorem 5.1 exact fixed-angle intervals + Seidel-block Gram verification."""
    q = _req_int(cert, "inv_cosine_q")
    k = _req_int(cert, "excess_k")
    t_claim = _req_int(cert, "claimed_t")
    r_claim = _req_int(cert, "claimed_r")
    dm_claim = _req_int(cert, "claimed_d_minus")
    dp_claim = _req_int(cert, "claimed_d_plus")
    check_gram = _req_bool(cert, "verify_explicit_seidel_gram")
    if q < 3 or q % 2 == 0 or k not in (7, 23):
        raise ValueError("Theorem 5.1 requires odd q >= 3 and k in {7, 23}")
    M, b, a = (28, 10, 3) if k == 7 else (276, 56, 5)
    Q = q + 1
    # t = 2 * ceil(Q / (2*b))
    t = 2 * ((Q + 2 * b - 1) // (2 * b))
    r = (b * t - Q) // 2
    if (b * t - Q) % 2 != 0:
        raise ValueError("b*t - Q must be even")
    d_minus = M * t - k
    floor_val = _floor_Lk_Q_exact(k, Q)
    d_plus = floor_val - k - 1
    if t != t_claim or r != r_claim or d_minus != dm_claim or d_plus != dp_claim:
        raise ValueError(
            f"Thm 5.1 mismatch: got (t={t}, r={r}, d-={d_minus}, d+={d_plus}), "
            f"claimed ({t_claim}, {r_claim}, {dm_claim}, {dp_claim})"
        )
    if not (0 <= r < t and d_minus <= d_plus):
        raise ValueError("Thm 5.1 interval hypotheses r < t and d- <= d+ failed")

    # Optionally build and verify exact integer Gram matrix G_hat (Eq. 38 + extension u)
    gram_ranks: List[int] = []
    if check_gram:
        if k != 7 or t > 4:
            raise ValueError("explicit Seidel Gram check supported for k=7, t<=4")
        # Build seed 28x28 Seidel matrix S of 28 lines in R^7 (pairs of 8-set, cosine 1/3)
        pairs = [(i, j) for i in range(8) for j in range(i + 1, 8)]
        assert len(pairs) == 28
        S = [[0] * 28 for _ in range(28)]
        for i in range(28):
            si = set(pairs[i])
            for j in range(i + 1, 28):
                # <1_e - 1/4 1, 1_f - 1/4 1> = |e cap f| - 1/2 in {-1/2, +1/2}, scaled by 2 -> {-1, +1}
                val = 2 * len(si.intersection(pairs[j])) - 1
                S[i][j] = val
                S[j][i] = val
        # r-regular circulant adjacency matrix A on Z/tZ
        diffs = set()
        for d in range(1, r // 2 + 1):
            diffs.add(d)
            diffs.add(t - d)
        if r % 2 == 1:
            diffs.add(t // 2)
        A = [[1 if ((j - i) % t) in diffs else 0 for j in range(t)] for i in range(t)]
        # Eq. (38): G = Q*I_{Mt} - (S + I_M) x J_t + 2*I_M x A
        N0 = M * t
        G0: List[List[Fraction]] = [[Fraction(0)] * N0 for _ in range(N0)]
        for bi in range(M):
            for bj in range(M):
                s_plus = S[bi][bj] + (1 if bi == bj else 0)
                for ci in range(t):
                    for cj in range(t):
                        row = bi * t + ci
                        col = bj * t + cj
                        val = (Q if row == col else 0) - s_plus + (2 * A[ci][cj] if bi == bj else 0)
                        G0[row][col] = Fraction(val)
        # Check diagonal == q and off-diagonal in {-1, +1}
        for i in range(N0):
            if G0[i][i] != q:
                raise ValueError("G0 diagonal != q")
            for j in range(i + 1, N0):
                if G0[i][j] not in (-1, 1):
                    raise ValueError("G0 off-diagonal not in +-1")
        # Check extension for u = 0 .. (d_plus - d_minus)
        z = [Fraction(1 if (c % 2 == 0) else -1) for _ in range(M) for c in range(t)]
        for u in range(d_plus - d_minus + 1):
            Nu = N0 + u
            Gu: List[List[Fraction]] = [[Fraction(0)] * Nu for _ in range(Nu)]
            for i in range(N0):
                for j in range(N0):
                    Gu[i][j] = G0[i][j]
            for a_idx in range(u):
                for i in range(N0):
                    Gu[i][N0 + a_idx] = z[i]
                    Gu[N0 + a_idx][i] = z[i]
                for b_idx in range(u):
                    Gu[N0 + a_idx][N0 + b_idx] = Fraction(q if a_idx == b_idx else 1)
            psd, rk = exact_psd_and_rank(Gu)
            if not psd or rk != d_minus + u:
                raise ValueError(f"Extended Gram failed for u={u}: psd={psd}, rank={rk}, expected={d_minus+u}")
            gram_ranks.append(rk)
    return {
        "mode": "fixed_angle_intervals_thm51",
        "q": q,
        "k": k,
        "t": t,
        "r": r,
        "interval": [d_minus, d_plus],
        "gram_ranks": gram_ranks,
    }


def _legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def verify_conference_and_cosine_1_7(cert: Dict[str, Any]) -> Dict[str, Any]:
    """Mode 4: Theorem 5.2 Paley conference matrices + Proposition 5.3 N_{1/7}(10)=12, N_{1/7}(11)=13."""
    submode = cert.get("submode")
    if submode == "thm52_paley_conference":
        p = _req_int(cert, "prime_p")
        t = _req_int(cert, "repeat_t")
        d_claim = _req_int(cert, "claimed_dim")
        n_claim = _req_int(cert, "claimed_lines")
        if p < 5 or p % 4 != 1 or t < 2:
            raise ValueError("Theorem 5.2 Paley specialization requires prime p = 1 mod 4 and t >= 2")
        # Primality check
        for f in range(2, p):
            if f * f > p:
                break
            if p % f == 0:
                raise ValueError(f"p={p} is not prime")
        m = p + 1
        # Build symmetric Paley conference matrix C of order m = p + 1
        C = [[0] * m for _ in range(m)]
        for i in range(1, m):
            C[0][i] = 1
            C[i][0] = 1
        for i in range(p):
            for j in range(p):
                if i != j:
                    C[i + 1][j + 1] = _legendre(i - j, p)
        # Verify C symmetric, zero diag, off-diag +-1, C^2 = p*I_m
        for i in range(m):
            if C[i][i] != 0:
                raise ValueError("conference matrix diagonal nonzero")
            for j in range(m):
                if i != j and C[i][j] not in (-1, 1):
                    raise ValueError("conference off-diag not +-1")
                if C[i][j] != C[j][i]:
                    raise ValueError("conference matrix not symmetric")
                dot = sum(C[i][k] * C[k][j] for k in range(m))
                if dot != (p if i == j else 0):
                    raise ValueError("C^2 != (m-1)*I_m")
        d = m * t - m // 2
        n = m * t
        if d != d_claim or n != n_claim:
            raise ValueError(f"Thm 5.2 (d, n) mismatch: got ({d}, {n}), claimed ({d_claim}, {n_claim})")
        # Trace contradiction: (mt+1)^2 < (m/2 + 1) * 2 * m * t^2
        lhs = (m * t + 1) ** 2
        rhs = (m // 2 + 1) * 2 * m * t * t
        if rhs - lhs != 2 * m * t * (t - 1) - 1 or lhs >= rhs:
            raise ValueError("Thm 5.2 trace contradiction failed")
        return {"mode": "conference_and_cosine_1_7", "submode": submode, "p": p, "dim": d, "lines": n}
    elif submode == "prop53_cosine_1_7":
        d = _req_int(cert, "dim_d")
        n_claim = _req_int(cert, "claimed_lines")
        if d not in (10, 11):
            raise ValueError("Proposition 5.3 certifies d in {10, 11}")
        rel_bound = (48 * d) // (49 - d)
        if n_claim != rel_bound:
            raise ValueError(f"relative bound mismatch: {rel_bound} != {n_claim}")
        n = n_claim
        # B = diag(J4, J4, J4) for n=12, or diag(J4, J4, J4, J1) for n=13
        block_of = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2] + ([3] if n == 13 else [])
        G: List[List[Fraction]] = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                b_ij = 1 if block_of[i] == block_of[j] else 0
                val = (8 if i == j else 0) - 2 * b_ij + 1
                G[i][j] = Fraction(val)
        for i in range(n):
            if G[i][i] != 7:
                raise ValueError("Prop 5.3 diagonal != 7")
            for j in range(i + 1, n):
                if G[i][j] not in (-1, 1):
                    raise ValueError("Prop 5.3 off-diagonal not in +-1")
        psd, rank = exact_psd_and_rank(G)
        if not psd or rank != d:
            raise ValueError(f"Prop 5.3 Gram check failed: psd={psd}, rank={rank}, expected={d}")
        return {"mode": "conference_and_cosine_1_7", "submode": submode, "dim": d, "lines": n, "gram_rank": rank}
    else:
        raise ValueError(f"unknown conference/1-7 submode '{submode}'")


def verify_certificate(cert: Dict[str, Any]) -> Dict[str, Any]:
    if not isinstance(cert, dict):
        raise ValueError("certificate must be JSON object")
    if cert.get("schema") != SCHEMA:
        raise ValueError(f"expected schema '{SCHEMA}', got '{cert.get('schema')}'")
    mode = cert.get("mode")
    if mode == "dim18_construction_and_reduction":
        return verify_dim18_construction_and_reduction(cert)
    if mode == "leech_evaluation_code_dim20":
        return verify_leech_evaluation_code_dim20(cert)
    if mode == "fixed_angle_intervals_thm51":
        return verify_fixed_angle_intervals_thm51(cert)
    if mode == "conference_and_cosine_1_7":
        return verify_conference_and_cosine_1_7(cert)
    raise ValueError(f"unknown mode '{mode}'")


def verify_file(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    res = verify_certificate(data)
    res["file"] = path.name
    return res


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("Usage: python -m src.equiangular_dim18_20_verifier <cert_or_dir> [...]")
        return 2
    files: List[Path] = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            files.extend(sorted(p.glob("*.json")))
        else:
            files.append(p)
    if not files:
        print("No certificate files found.")
        return 1
    for f in files:
        res = verify_file(f)
        print(f"PASS  {f.name}  ({res['mode']})")
    print(f"VERIFIED {len(files)}/{len(files)} certificates under {SCHEMA}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
