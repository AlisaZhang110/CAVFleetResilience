"""
Stage I -- PTN fault isolation.

Implements Eqs. (1)-(25) of Liu et al. (2025).

Equation-by-equation map (deviations flagged DEV-n, explained in README):

  (1)(2)  d_i >= mu_ij * a_ij,0                 endpoints of a faulted line are faulted
  (3)(4)  d_i + a_ij,I - 1 <= d_j               zone spreads through *still-closed* lines
                                                DEV-1: paper prints a_ij,0 here; with the
                                                printed form no switching action could ever
                                                contain a zone, and Case 1 vs Case 2 would
                                                be identical.  Stage II's Eqs. (29)-(30) use
                                                the post-action status, so a_I is intended.
  (5)     d_utn_node >= d_pdn_bus  at each FCS  intra-propagation, PDN -> UTN only
                                                DEV-2: paper prints an equality, which also
                                                forces UTN -> PDN; its own text says that
                                                direction is excluded.  Implemented as the
                                                one-way implication the text describes.
  (6)     faulted UTN link capacity -> 0
  (7)(8)  folded into (9): a switch lets a line be opened
  (9)     (1-mu)(a0 - k) <= a_I <= (1-mu) a0    Stage I may only OPEN lines, never close
  (10)-(20) linearised DistFlow with line status, load shedding, generator limits
  (21)    a_ij,I = 1 - e_i(j)                   DEV-3: paper prints a = e, which inverts the
                                                meaning (a=1 is "unobstructed", e=1 is
                                                "blocked").  Sign corrected.
  (22)    sum_j e_i(j) <= N_i - 1               DEV-4: paper prints <= N_c (undefined symbol,
                                                and vacuous).  A budget of N_i - 1 keeps at
                                                least one approach open at every crossing.
  (23)(24) blocked link capacity -> 0
  (25)    objective                             DEV-5: the two printed terms have
                                                incompatible units (MWh-weighted load vs.
                                                vehicle capacity) and are summed unweighted.
                                                Implemented with the UTN term normalised to
                                                a failure fraction and scaled by LAMBDA_CAP.

Usable UTN capacity additionally requires both endpoints to be non-faulted
(Q* <= Q0 (1 - d_i), Q* <= Q0 (1 - d_j)).  This is not a numbered equation in
the paper but it is what makes the Case 1 "failure rate = 100%" statement
arithmetically possible: without it a fully-propagated zone would still report
full capacity.
"""
from __future__ import annotations

from . import data as D
from .lp import Model, INF


class Stage1Result:
    def __init__(self):
        self.d_pdn = {}
        self.d_utn = {}
        self.a_pdn = {}
        self.a_utn = {}
        self.qstar = {}          # directed UTN link -> capacity after Stage I
        self.pL = {}
        self.load_loss_weighted = 0.0
        self.failure_rate = 0.0
        self.obj = 0.0
        self.stats = {}


def solve_stage1(tb, allow_crossing_elimination: bool, verbose=False) -> Stage1Result:
    m = Model("stage1")
    BIGM_V = 2.0

    # ---------------- variables ----------------
    d_p = {b: m.binary(f"dP{b}") for b in tb.pdn_buses}
    d_t = {n: m.binary(f"dT{n}") for n in tb.utn_nodes}
    a_p = {ln: m.binary(f"aP{ln}") for ln in tb.pdn_lines}
    a_t = {ln: m.binary(f"aT{ln}") for ln in tb.utn_links}
    e_t = {ln: m.binary(f"e{ln}") for ln in tb.utn_links}     # block link (i,j) at crossing i
    qs = {ln: m.var(f"Q*{ln}", 0.0, tb.q0[ln]) for ln in tb.utn_links}

    P = {ln: m.free(f"P{ln}") for ln in tb.pdn_lines}
    Q = {ln: m.free(f"Q{ln}") for ln in tb.pdn_lines}
    u = {b: m.var(f"u{b}", D.VMIN, D.VMAX) for b in tb.pdn_buses}
    pL = {b: m.var(f"pL{b}", 0.0, max(tb.pl[b], 1e-9)) for b in tb.pdn_buses}
    qL = {b: m.var(f"qL{b}", 0.0, max(tb.ql[b], 1e-9)) for b in tb.pdn_buses}
    pg = {b: m.var(f"pg{b}", 0.0, tb.gen.get(b, (0, 0))[0]) for b in tb.pdn_buses}
    qg = {b: m.var(f"qg{b}", 0.0, tb.gen.get(b, (0, 0))[1]) for b in tb.pdn_buses}

    m.constr({u[D.SUBSTATION_BUS]: 1.0}, "==", D.VREF, "slack voltage")

    # ---------------- (9) line status ----------------
    for ln in tb.pdn_lines:
        mu, a0, k = tb.mu_pdn[ln], tb.a0_pdn[ln], tb.k0_pdn[ln]
        m.constr({a_p[ln]: 1.0}, "<=", (1 - mu) * a0, f"eq9u{ln}")
        m.constr({a_p[ln]: 1.0}, ">=", (1 - mu) * (a0 - k), f"eq9l{ln}")

    for ln in tb.utn_links:
        mu = tb.mu_utn[ln]
        # (21) corrected sign, folded with the fault status
        m.constr({a_t[ln]: 1.0, e_t[ln]: 1.0}, "==", 1 - mu, f"eq21{ln}")
        if not allow_crossing_elimination:
            m.constr({e_t[ln]: 1.0}, "==", 0.0, f"nocross{ln}")

    # (22) budget of roadblocks per crossing.
    # A healthy crossing must keep at least one approach open; a faulted
    # crossing may be sealed completely, which is exactly the action the
    # paper describes ("roadblocks at the critical intersections" so that
    # vehicles cannot enter links connected to faulted nodes).
    if allow_crossing_elimination:
        by_node = {}
        for (i, j) in tb.utn_links:
            by_node.setdefault(i, []).append((i, j))
        for i, lst in by_node.items():
            row = {e_t[ln]: 1.0 for ln in lst}
            row[d_t[i]] = -float(len(lst))
            m.constr(row, "<=", len(lst) - 1, f"eq22{i}")

    # ---------------- (1)-(4) inter-propagation ----------------
    for ln in tb.pdn_lines:
        i, j = ln
        mu, a0 = tb.mu_pdn[ln], tb.a0_pdn[ln]
        m.constr({d_p[i]: 1.0}, ">=", mu * a0, f"eq1{ln}")
        m.constr({d_p[j]: 1.0}, ">=", mu * a0, f"eq2{ln}")
        m.constr({d_p[i]: 1.0, a_p[ln]: 1.0, d_p[j]: -1.0}, "<=", 1.0, f"eq3{ln}")
        m.constr({d_p[j]: 1.0, a_p[ln]: 1.0, d_p[i]: -1.0}, "<=", 1.0, f"eq4{ln}")

    for ln in tb.utn_links:
        i, j = ln
        mu = tb.mu_utn[ln]
        m.constr({d_t[i]: 1.0}, ">=", float(mu), f"eq1t{ln}")
        m.constr({d_t[j]: 1.0}, ">=", float(mu), f"eq2t{ln}")
        m.constr({d_t[i]: 1.0, a_t[ln]: 1.0, d_t[j]: -1.0}, "<=", 1.0, f"eq3t{ln}")
        m.constr({d_t[j]: 1.0, a_t[ln]: 1.0, d_t[i]: -1.0}, "<=", 1.0, f"eq4t{ln}")

    # ---------------- (5) intra-propagation PDN -> UTN ----------------
    for f, bus in tb.fcs_bus.items():
        node = tb.fcs_node[f]
        m.constr({d_t[node]: 1.0, d_p[bus]: -1.0}, ">=", 0.0, f"eq5fcs{f}")

    # ---------------- (6)(24) UTN capacity ----------------
    for ln in tb.utn_links:
        i, j = ln
        q0 = tb.q0[ln]
        m.constr({qs[ln]: 1.0, a_t[ln]: -q0}, "<=", 0.0, f"eq24{ln}")
        m.constr({qs[ln]: 1.0, d_t[i]: q0}, "<=", q0, f"cap_di{ln}")
        m.constr({qs[ln]: 1.0, d_t[j]: q0}, "<=", q0, f"cap_dj{ln}")

    # ---------------- (10)-(20) PDN operation ----------------
    inc = {b: [] for b in tb.pdn_buses}
    out = {b: [] for b in tb.pdn_buses}
    for ln in tb.pdn_lines:
        i, j = ln
        out[i].append(ln)
        inc[j].append(ln)

    for b in tb.pdn_buses:
        row = {pL[b]: -1.0, pg[b]: -1.0}
        for ln in inc[b]:
            row[P[ln]] = row.get(P[ln], 0.0) - 1.0
        for ln in out[b]:
            row[P[ln]] = row.get(P[ln], 0.0) + 1.0
        m.constr(row, "==", -tb.pl[b], f"eq10bus{b}")

        row = {qL[b]: -1.0, qg[b]: -1.0}
        for ln in inc[b]:
            row[Q[ln]] = row.get(Q[ln], 0.0) - 1.0
        for ln in out[b]:
            row[Q[ln]] = row.get(Q[ln], 0.0) + 1.0
        m.constr(row, "==", -tb.ql[b], f"eq11bus{b}")

    for ln in tb.pdn_lines:
        i, j = ln
        rr, xx = tb.r[ln] / D.SBASE, tb.x[ln] / D.SBASE
        m.constr({u[i]: 1.0, u[j]: -1.0, P[ln]: -rr, Q[ln]: -xx, a_p[ln]: BIGM_V},
                 "<=", BIGM_V, f"eq12{ln}")
        m.constr({u[i]: 1.0, u[j]: -1.0, P[ln]: -rr, Q[ln]: -xx, a_p[ln]: -BIGM_V},
                 ">=", -BIGM_V, f"eq13{ln}")
        s = tb.smax[ln]
        m.range_constr({P[ln]: 1.0, a_p[ln]: -s}, -INF, 0.0, f"eq14u{ln}")
        m.range_constr({P[ln]: 1.0, a_p[ln]: s}, 0.0, INF, f"eq14l{ln}")
        m.range_constr({Q[ln]: 1.0, a_p[ln]: -s}, -INF, 0.0, f"eq15u{ln}")
        m.range_constr({Q[ln]: 1.0, a_p[ln]: s}, 0.0, INF, f"eq15l{ln}")

    for b, (pmax, qmax) in tb.gen.items():
        m.constr({pg[b]: 1.0, d_p[b]: pmax}, "<=", pmax, f"eq17{b}")
        m.constr({qg[b]: 1.0, d_p[b]: qmax}, "<=", qmax, f"eq18{b}")

    for b in tb.pdn_buses:
        if tb.pl[b] > 0:
            m.constr({pL[b]: 1.0, d_p[b]: -tb.pl[b]}, ">=", 0.0, f"eq19{b}")
        if tb.ql[b] > 0:
            m.constr({qL[b]: 1.0, d_p[b]: -tb.ql[b]}, ">=", 0.0, f"eq20{b}")

    # ---------------- (25) objective ----------------
    dur = D.TI_H - D.T0_H
    obj = {pL[b]: dur * tb.w[b] for b in tb.pdn_buses}
    qtot = sum(D.ZETA_LINK * tb.q0[ln] for ln in tb.utn_links)
    for ln in tb.utn_links:
        obj[qs[ln]] = obj.get(qs[ln], 0.0) - D.LAMBDA_CAP_STAGE1 * D.ZETA_LINK / qtot
    m.minimize(obj)

    sol = m.solve(verbose=verbose)

    res = Stage1Result()
    res.obj = sol.obj
    res.stats = dict(n_vars=m.n_vars, n_int=m.n_int, n_constr=m.n_constr)
    res.d_pdn = {b: sol.b(v) for b, v in d_p.items()}
    res.d_utn = {n: sol.b(v) for n, v in d_t.items()}
    res.a_pdn = {ln: sol.b(v) for ln, v in a_p.items()}
    res.a_utn = {ln: sol.b(v) for ln, v in a_t.items()}
    res.qstar = {ln: max(0.0, sol[v]) for ln, v in qs.items()}
    res.pL = {b: sol[v] for b, v in pL.items()}
    res.load_loss_weighted = sum(tb.w[b] * res.pL[b] for b in tb.pdn_buses)
    res.failure_rate = 1.0 - sum(D.ZETA_LINK * res.qstar[ln] for ln in tb.utn_links) / qtot
    return res
