"""
Stage II -- PTN service restoration.

Implements Eqs. (26)-(60) of Liu et al. (2025).

Equation map (deviations DEV-n explained in README):

  (26)      (1-mu)(a_I - k) <= a_II <= (1-mu)(a_I + k)   Stage II may close tie lines
  (27)(28)  z_ij + z_ji = a_II ; sum_in z <= 1 - k_j - dk_j gamma_j   radiality
  (29)(30)  no reconnection across a faulted-zone boundary
  (31)-(41) multi-period linearised DistFlow, shedding, generator limits
  (42)(43)  coupled FCS bus:  p_load = p_FCS + p_0
  (44)-(46) link contraflow.  Q*_I is a *parameter* here (Stage I output), so
            Eq. (44) is linear as printed.
  (47)(48)  expanded transportation graph E = (N, L), N subset of N_T x T x C
            DEV-6: the printed L_R forces a vehicle to traverse a link every
            period (no dwell option).  A zero-cost waiting link is added,
            without which the fleet model is infeasible for any horizon longer
            than the shortest path.
  (49)      BPR-style congestion.  DEV-7: as printed t* = beta t0 (f/Q), i.e.
            the classical BPR exponent and the +1 term are dropped, and beta is
            simultaneously "economic value of time".  Kept as printed, but note
            that f * t*(f) is then QUADRATIC in f, so the printed Stage II
            problem is an MIQP, not the MILP the text implies.  Handled here by
            an outer (tangent-plane) linearisation of the convex quadratic.
  (50)-(53) customer trip requests, occupied AEV departures/arrivals
  (54)      link capacity.  DEV-8: as printed only *unoccupied* flow f_e is
            counted against capacity and cost, while occupied AEVs are invisible
            to congestion.  That cannot be the intent -- the whole contraflow
            result is about speeding up an occupied trip -- so occupied flows are
            counted here too.
  (55)      FCS capacity
  (57)      p_FCS = f_e * (c_n - c_m);  >0 charge (load), <0 discharge (V2G)
  (58)(59)(60) system-optimal objective with weights theta_P, theta_T

Unserved trip demand is priced at a large penalty rather than made infeasible,
because Case 1 (100% UTN failure) genuinely admits no feasible routing and the
paper's own Table III silently omits it.
"""
from __future__ import annotations

from . import data as D
from .lp import Model, INF

UNSERVED_PENALTY = 1.0e4


class Stage2Result:
    def __init__(self):
        self.a_pdn = {}
        self.contraflow = {}
        self.qstar = {}
        self.pL = {}            # (bus, t) -> shed MW
        self.pfcs = {}          # (fcs, t) -> MW  (+ charge, - discharge)
        self.travel_cost = 0.0
        self.load_loss_weighted = 0.0   # MWh-weighted
        self.unserved = 0.0
        self.obj = 0.0
        self.stats = {}


def _expanded_graph(tb, s1, allow_v2g, faulted_nodes):
    """Build E = (N, L).  Returns (links, meta) where meta[k] describes link k."""
    T = D.N_TIME_NODES
    C = D.SOC_LEVELS
    links = []          # (tail_node, head_node)
    meta = []           # dict(kind=..., utn=..., t=..., dsoc=...)

    def node(n, t, c):
        return (n, t, c)

    for t in range(1, T):
        # --- travel links (L_R) ---
        for (i, j) in tb.utn_links:
            if i in faulted_nodes or j in faulted_nodes:
                continue
            for c in range(1 + D.TRAVEL_SOC_COST, C + 1):
                links.append((node(i, t, c), node(j, t + 1, c - D.TRAVEL_SOC_COST)))
                meta.append(dict(kind="R", utn=(i, j), t=t, dsoc=-D.TRAVEL_SOC_COST))
        # --- waiting links (DEV-6) ---
        for n in tb.utn_nodes:
            if n in faulted_nodes:
                continue
            for c in range(1, C + 1):
                links.append((node(n, t, c), node(n, t + 1, c)))
                meta.append(dict(kind="W", utn=None, t=t, dsoc=0))
        # --- FCS charging / discharging links (L_FCS) ---
        for f, nd in tb.fcs_node.items():
            if nd in faulted_nodes or tb.fcs_bus[f] in s1.d_pdn and s1.d_pdn[tb.fcs_bus[f]]:
                continue
            for c in range(1, C):
                links.append((node(nd, t, c), node(nd, t + 1, c + 1)))
                meta.append(dict(kind="FCS", utn=None, t=t, dsoc=+1, fcs=f))
            if allow_v2g:
                for c in range(2, C + 1):
                    links.append((node(nd, t, c), node(nd, t + 1, c - 1)))
                    meta.append(dict(kind="FCS", utn=None, t=t, dsoc=-1, fcs=f))
    return links, meta


def solve_stage2(tb, s1, allow_contraflow, allow_v2g, theta_t_over_p=1.0,
                 trip_requests=None, verbose=False, n_tangents=6) -> Stage2Result:
    m = Model("stage2")
    T = D.N_TIME_NODES
    periods = list(range(1, D.N_PERIODS + 1))
    dt = (D.TII_H - D.TI_H) / D.N_PERIODS
    trips = trip_requests if trip_requests is not None else D.TRIP_REQUESTS
    BIGM_V = 2.0

    faulted_nodes = {n for n, v in s1.d_utn.items() if v}
    # links that Stage I left blocked are unusable in Stage II as well
    blocked = {ln for ln, v in s1.a_utn.items() if not v}

    # ================= PDN =================
    a = {ln: m.binary(f"a2{ln}") for ln in tb.pdn_lines}
    z = {}
    for ln in tb.pdn_lines:
        i, j = ln
        z[(i, j)] = m.binary(f"z{i}_{j}")
        z[(j, i)] = m.binary(f"z{j}_{i}")
    gam = {b: m.binary(f"g{b}") for b in D.DG_BUSES}

    for ln in tb.pdn_lines:
        mu, k = tb.mu_pdn[ln], tb.k0_pdn[ln]
        aI = 1 if s1.a_pdn[ln] else 0
        m.constr({a[ln]: 1.0}, "<=", (1 - mu) * (aI + k), f"eq26u{ln}")
        m.constr({a[ln]: 1.0}, ">=", (1 - mu) * (aI - k), f"eq26l{ln}")
        # (29)(30): cannot bridge a faulted-zone boundary
        di, dj = int(s1.d_pdn[ln[0]]), int(s1.d_pdn[ln[1]])
        if di != dj:
            m.constr({a[ln]: 1.0}, "==", 0.0, f"eq2930{ln}")
        # (27)
        m.constr({z[(ln[0], ln[1])]: 1.0, z[(ln[1], ln[0])]: 1.0, a[ln]: -1.0},
                 "==", 0.0, f"eq27{ln}")

    inc_z = {b: [] for b in tb.pdn_buses}
    for ln in tb.pdn_lines:
        i, j = ln
        inc_z[j].append(z[(i, j)])
        inc_z[i].append(z[(j, i)])
    for b in tb.pdn_buses:
        row = {v: 1.0 for v in inc_z[b]}
        kb = 1.0 if b == D.SUBSTATION_BUS else 0.0
        rhs = 1.0 - kb
        if b in D.DG_BUSES:
            row[gam[b]] = 1.0
        m.constr(row, "<=", rhs, f"eq28{b}")

    P, Qf, u, pL, qL, pg, qg = {}, {}, {}, {}, {}, {}, {}
    for t in periods:
        for ln in tb.pdn_lines:
            P[(ln, t)] = m.free(f"P{ln}_{t}")
            Qf[(ln, t)] = m.free(f"Q{ln}_{t}")
        for b in tb.pdn_buses:
            u[(b, t)] = m.var(f"u{b}_{t}", D.VMIN, D.VMAX)
            pL[(b, t)] = m.var(f"pL{b}_{t}", 0.0, max(tb.pl[b], 1e-9))
            qL[(b, t)] = m.var(f"qL{b}_{t}", 0.0, max(tb.ql[b], 1e-9))
            pg[(b, t)] = m.var(f"pg{b}_{t}", 0.0, tb.gen.get(b, (0, 0))[0])
            qg[(b, t)] = m.var(f"qg{b}_{t}", 0.0, tb.gen.get(b, (0, 0))[1])
        m.constr({u[(D.SUBSTATION_BUS, t)]: 1.0}, "==", D.VREF, f"slack{t}")

    # FCS net power (42)(43)(57): free sign, fixed later by AEV flows
    pfcs = {(f, t): m.free(f"pFCS{f}_{t}") for f in tb.fcs_bus for t in periods}

    inc = {b: [] for b in tb.pdn_buses}
    out = {b: [] for b in tb.pdn_buses}
    for ln in tb.pdn_lines:
        out[ln[0]].append(ln)
        inc[ln[1]].append(ln)
    bus_fcs = {}
    for f, b in tb.fcs_bus.items():
        bus_fcs.setdefault(b, []).append(f)

    for t in periods:
        for b in tb.pdn_buses:
            row = {pL[(b, t)]: -1.0, pg[(b, t)]: -1.0}
            for ln in inc[b]:
                row[P[(ln, t)]] = row.get(P[(ln, t)], 0.0) - 1.0
            for ln in out[b]:
                row[P[(ln, t)]] = row.get(P[(ln, t)], 0.0) + 1.0
            for f in bus_fcs.get(b, []):
                row[pfcs[(f, t)]] = row.get(pfcs[(f, t)], 0.0) + 1.0
            m.constr(row, "==", -tb.pl[b], f"eq31_{b}_{t}")

            row = {qL[(b, t)]: -1.0, qg[(b, t)]: -1.0}
            for ln in inc[b]:
                row[Qf[(ln, t)]] = row.get(Qf[(ln, t)], 0.0) - 1.0
            for ln in out[b]:
                row[Qf[(ln, t)]] = row.get(Qf[(ln, t)], 0.0) + 1.0
            m.constr(row, "==", -tb.ql[b], f"eq32_{b}_{t}")

        for ln in tb.pdn_lines:
            i, j = ln
            rr, xx = tb.r[ln] / D.SBASE, tb.x[ln] / D.SBASE
            m.constr({u[(i, t)]: 1.0, u[(j, t)]: -1.0, P[(ln, t)]: -rr,
                      Qf[(ln, t)]: -xx, a[ln]: BIGM_V}, "<=", BIGM_V, f"eq33_{ln}_{t}")
            m.constr({u[(i, t)]: 1.0, u[(j, t)]: -1.0, P[(ln, t)]: -rr,
                      Qf[(ln, t)]: -xx, a[ln]: -BIGM_V}, ">=", -BIGM_V, f"eq34_{ln}_{t}")
            s = tb.smax[ln]
            m.range_constr({P[(ln, t)]: 1.0, a[ln]: -s}, -INF, 0.0, f"eq35u{ln}{t}")
            m.range_constr({P[(ln, t)]: 1.0, a[ln]: s}, 0.0, INF, f"eq35l{ln}{t}")
            m.range_constr({Qf[(ln, t)]: 1.0, a[ln]: -s}, -INF, 0.0, f"eq36u{ln}{t}")
            m.range_constr({Qf[(ln, t)]: 1.0, a[ln]: s}, 0.0, INF, f"eq36l{ln}{t}")

        for b, (pmax, qmax) in tb.gen.items():
            if s1.d_pdn[b]:
                m.constr({pg[(b, t)]: 1.0}, "==", 0.0, f"eq38{b}{t}")
                m.constr({qg[(b, t)]: 1.0}, "==", 0.0, f"eq39{b}{t}")

        for b in tb.pdn_buses:
            if s1.d_pdn[b]:
                if tb.pl[b] > 0:
                    m.constr({pL[(b, t)]: 1.0}, "==", tb.pl[b], f"eq40{b}{t}")
                if tb.ql[b] > 0:
                    m.constr({qL[(b, t)]: 1.0}, "==", tb.ql[b], f"eq41{b}{t}")

    # ================= UTN contraflow (44)-(46) =================
    cflow = {ln: m.binary(f"c{ln}") for ln in tb.utn_links}
    qs2 = {}
    for ln in tb.utn_links:
        i, j = ln
        rev = (j, i)
        qI, qIrev = s1.qstar[ln], s1.qstar[rev]
        if not allow_contraflow:
            m.constr({cflow[ln]: 1.0}, "==", 0.0, f"nocf{ln}")
        v = m.var(f"Q2{ln}", 0.0, qI + qIrev + 1e-9)
        qs2[ln] = v
        # Q*_ij,II = Q*_ij,I (1 - c_ij) + Q*_ji,I c_ji
        m.constr({v: 1.0, cflow[ln]: qI, cflow[rev]: -qIrev}, "==", qI, f"eq44{ln}")
    for (i, j) in D.utn_undirected(tb):
        m.constr({cflow[(i, j)]: 1.0, cflow[(j, i)]: 1.0}, "<=", 1.0, f"eq45{i}_{j}")
    m.constr({v: 1.0 for v in cflow.values()}, "<=", float(len(tb.utn_links)), "eq46")

    # ================= AEV fleet =================
    links, meta = _expanded_graph(tb, s1, allow_v2g, faulted_nodes)
    keep = [k for k in range(len(links))
            if meta[k]["utn"] is None or (meta[k]["utn"] not in blocked)]
    links = [links[k] for k in keep]
    meta = [meta[k] for k in keep]

    fu = [m.var(f"fu{k}") for k in range(len(links))]
    fo = {}
    for oi, _ in enumerate(trips):
        fo[oi] = [m.var(f"fo{oi}_{k}") if meta[k]["kind"] in ("R", "W") else None
                  for k in range(len(links))]

    nodes = set()
    for (tl, hd) in links:
        nodes.add(tl)
        nodes.add(hd)
    outk = {v: [] for v in nodes}
    ink = {v: [] for v in nodes}
    for k, (tl, hd) in enumerate(links):
        outk[tl].append(k)
        ink[hd].append(k)

    # initial / terminal fleet
    # Vehicles initially parked at a node that Stage I sealed off are stranded
    # and drop out of the Stage II fleet; count only what the graph can hold.
    Sf = {}
    for (n, c, cnt) in D.FLEET_INIT:
        if (n, 1, c) in nodes:
            Sf[(n, 1, c)] = Sf.get((n, 1, c), 0.0) + cnt
    fleet_total = sum(Sf.values())
    Ef = {v: m.var(f"Ef{v}") for v in nodes if v[1] == T}

    # occupied departures / arrivals (50)(52)(53)
    dep, arr, unserved = {}, {}, {}
    for oi, (mo, no, to, alpha) in enumerate(trips):
        unserved[oi] = m.var(f"uns{oi}", 0.0, alpha)
        dep[oi] = {}
        arr[oi] = {}
        for c in range(1, D.SOC_LEVELS + 1):
            if (mo, to, c) in nodes:
                dep[oi][c] = m.var(f"dep{oi}_{c}")
        for t in range(to + 1, T + 1):
            for c in range(1, D.SOC_LEVELS + 1):
                if (no, t, c) in nodes:
                    arr[oi][(t, c)] = m.var(f"arr{oi}_{t}_{c}")
        m.constr({**{v: 1.0 for v in dep[oi].values()}, unserved[oi]: 1.0},
                 "==", alpha, f"eq52_{oi}")
        m.constr({**{v: 1.0 for v in arr[oi].values()}, unserved[oi]: 1.0},
                 "==", alpha, f"eq53_{oi}")

    # commodity balance for occupied flows
    for oi, (mo, no, to, alpha) in enumerate(trips):
        for v in nodes:
            row = {}
            for k in ink[v]:
                if fo[oi][k] is not None:
                    row[fo[oi][k]] = row.get(fo[oi][k], 0.0) + 1.0
            for k in outk[v]:
                if fo[oi][k] is not None:
                    row[fo[oi][k]] = row.get(fo[oi][k], 0.0) - 1.0
            if v in dep[oi] or (v[0] == mo and v[1] == to and v[2] in dep[oi]):
                pass
            if v[0] == mo and v[1] == to and v[2] in dep[oi]:
                row[dep[oi][v[2]]] = row.get(dep[oi][v[2]], 0.0) + 1.0
            if v[0] == no and (v[1], v[2]) in arr[oi]:
                row[arr[oi][(v[1], v[2])]] = row.get(arr[oi][(v[1], v[2])], 0.0) - 1.0
            m.constr(row, "==", 0.0, f"occ{oi}_{v}")

    # unoccupied balance (51)
    for v in nodes:
        row = {}
        for k in ink[v]:
            row[fu[k]] = row.get(fu[k], 0.0) + 1.0
        for k in outk[v]:
            row[fu[k]] = row.get(fu[k], 0.0) - 1.0
        for oi, (mo, no, to, alpha) in enumerate(trips):
            if v[0] == mo and v[1] == to and v[2] in dep[oi]:
                row[dep[oi][v[2]]] = row.get(dep[oi][v[2]], 0.0) - 1.0
            if v[0] == no and (v[1], v[2]) in arr[oi]:
                row[arr[oi][(v[1], v[2])]] = row.get(arr[oi][(v[1], v[2])], 0.0) + 1.0
        if v in Ef:
            row[Ef[v]] = row.get(Ef[v], 0.0) - 1.0
        m.constr(row, "==", -Sf.get(v, 0.0), f"eq51_{v}")
    m.constr({v: 1.0 for v in Ef.values()}, "==", float(fleet_total), "fleet_end")

    # link capacity (54) and per-(link,period) total flow
    Ftot = {}
    for ln in tb.utn_links:
        for t in periods:
            ks = [k for k in range(len(links))
                  if meta[k]["kind"] == "R" and meta[k]["utn"] == ln and meta[k]["t"] == t]
            if not ks:
                continue
            F = m.var(f"F{ln}_{t}", 0.0, D.Q0_LINK * 2)
            Ftot[(ln, t)] = F
            row = {F: -1.0}
            for k in ks:
                row[fu[k]] = row.get(fu[k], 0.0) + 1.0
                for oi in fo:
                    if fo[oi][k] is not None:
                        row[fo[oi][k]] = row.get(fo[oi][k], 0.0) + 1.0
            m.constr(row, "==", 0.0, f"Fdef{ln}_{t}")
            m.constr({F: 1.0, qs2[ln]: -1.0}, "<=", 0.0, f"eq54{ln}_{t}")

    # FCS capacity (55)
    for f, nd in tb.fcs_node.items():
        for t in periods:
            ks = [k for k in range(len(links))
                  if meta[k]["kind"] == "FCS" and meta[k].get("fcs") == f and meta[k]["t"] == t]
            if ks:
                m.constr({fu[k]: 1.0 for k in ks}, "<=", D.FCS_CAP, f"eq55{f}_{t}")

    # (57) FCS power coupling
    for f in tb.fcs_bus:
        for t in periods:
            ks = [k for k in range(len(links))
                  if meta[k]["kind"] == "FCS" and meta[k].get("fcs") == f and meta[k]["t"] == t]
            row = {pfcs[(f, t)]: 1.0}
            for k in ks:
                row[fu[k]] = row.get(fu[k], 0.0) - meta[k]["dsoc"] * D.SOC_STEP_MW
            m.constr(row, "==", 0.0, f"eq57{f}_{t}")

    # (49)(58) congestion cost, tangent-plane outer linearisation of beta*t0*F^2/Q
    zc = {}
    BIGM_C = 1.0e4
    for ln in tb.utn_links:
        i, j = ln
        rev = (j, i)
        qI, qIrev = s1.qstar[ln], s1.qstar[rev]
        for t in periods:
            if (ln, t) not in Ftot:
                continue
            F = Ftot[(ln, t)]
            v = m.var(f"zc{ln}_{t}", 0.0, BIGM_C)
            zc[(ln, t)] = v
            for cap, gate in ((qI, ("off", cflow[rev])), (qI + qIrev, ("on", cflow[rev]))):
                if cap <= 1e-9:
                    continue
                coef = D.BETA_VOT * tb.t0[ln] / cap
                for s_ in range(n_tangents + 1):
                    Fk = cap * s_ / n_tangents
                    # z >= coef (2 Fk F - Fk^2)  - M * (gate is inactive)
                    row = {v: 1.0, F: -2.0 * coef * Fk}
                    if gate[0] == "off":
                        row[gate[1]] = BIGM_C
                        rhs = -coef * Fk * Fk
                    else:
                        row[gate[1]] = -BIGM_C
                        rhs = -coef * Fk * Fk - BIGM_C
                    m.constr(row, ">=", rhs, f"bpr{ln}_{t}_{gate[0]}_{s_}")

    # ================= objective (60) =================
    theta_p = 1.0
    theta_t = theta_t_over_p
    obj = {}
    for t in periods:
        for b in tb.pdn_buses:
            obj[pL[(b, t)]] = obj.get(pL[(b, t)], 0.0) + theta_p * dt * tb.w[b]
    for v in zc.values():
        obj[v] = obj.get(v, 0.0) + theta_t
    for oi in unserved:
        # must dominate any achievable travel cost, including large theta_T
        obj[unserved[oi]] = UNSERVED_PENALTY * max(1.0, theta_t)
    m.minimize(obj)

    sol = m.solve(verbose=verbose, time_limit=600.0)

    res = Stage2Result()
    res.obj = sol.obj
    res.stats = dict(n_vars=m.n_vars, n_int=m.n_int, n_constr=m.n_constr,
                     n_expanded_links=len(links))
    res.a_pdn = {ln: sol.b(v) for ln, v in a.items()}
    res.contraflow = {ln: sol.b(v) for ln, v in cflow.items() if sol.b(v)}
    res.qstar = {ln: sol[v] for ln, v in qs2.items()}
    res.pL = {(b, t): sol[pL[(b, t)]] for b in tb.pdn_buses for t in periods}
    res.pfcs = {(f, t): sol[pfcs[(f, t)]] for f in tb.fcs_bus for t in periods}
    # Evaluate Eq. (49)/(58) directly on the realised flows rather than reading
    # the epigraph variable: with theta_T = 0 the epigraph is unpriced and would
    # float to its upper bound.
    tc = 0.0
    for (ln, t), Fv in Ftot.items():
        F = sol[Fv]
        cap = sol[qs2[ln]]
        if F > 1e-9 and cap > 1e-9:
            tc += D.BETA_VOT * tb.t0[ln] * F * F / cap
    res.travel_cost = tc
    res.travel_cost_epigraph = sum(sol[v] for v in zc.values())
    res.flows = {(ln, t): sol[v] for (ln, t), v in Ftot.items() if sol[v] > 1e-6}
    res.load_loss_weighted = sum(dt * tb.w[b] * res.pL[(b, t)]
                                 for b in tb.pdn_buses for t in periods)
    res.unserved = sum(sol[v] for v in unserved.values())
    return res
