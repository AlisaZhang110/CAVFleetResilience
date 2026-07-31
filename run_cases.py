"""
Reproduce the case study of Liu et al. (2025), Sec. IV, on the P33T12 testbed.

    python3 run_cases.py            # Cases 1-5 + resilience metrics
    python3 run_cases.py --sweep    # also the theta_T/theta_P sensitivity (Table VI)
    python3 run_cases.py --plots    # also write figures to out/

Case definitions follow Table II:

    Case   PDN recfg S1  PDN recfg S2  UTN recfg S1  UTN recfg S2  V2G  RO
      1         y             y             -             -         -    -
      2         y             y             y             -         -    -
      3         y             y             y             y         -    -
      4         y             y             y             y         y    -
      5         y             y             y             y         y    y
"""
from __future__ import annotations

import sys

from ptn import data as D
from ptn.stage1 import solve_stage1
from ptn.stage2 import solve_stage2

CASES = {
    1: dict(cross=False, contra=False, v2g=False, robust=False),
    2: dict(cross=True,  contra=False, v2g=False, robust=False),
    3: dict(cross=True,  contra=True,  v2g=False, robust=False),
    4: dict(cross=True,  contra=True,  v2g=True,  robust=False),
    5: dict(cross=True,  contra=True,  v2g=True,  robust=True),
}


def resilience_metric(tb, s1, s2):
    """Eq. (62).  Re = 1 - (weighted energy not served) / (weighted energy demanded)."""
    dur1 = D.TI_H - D.T0_H
    dt = (D.TII_H - D.TI_H) / D.N_PERIODS
    num = dur1 * sum(tb.w[b] * s1.pL[b] for b in tb.pdn_buses)
    num += dt * sum(tb.w[b] * s2.pL[(b, t)]
                    for b in tb.pdn_buses for t in range(1, D.N_PERIODS + 1))
    den = (D.TII_H - D.T0_H) * sum(tb.w[b] * tb.pl[b] for b in tb.pdn_buses)
    return 1.0 - num / den


def run_case(tb, cid, theta=None):
    cfg = CASES[cid]
    s1 = solve_stage1(tb, allow_crossing_elimination=cfg["cross"])
    trips = D.TRIP_REQUESTS
    if cfg["robust"]:
        trips = [(o, d, t, a * D.ROBUST_ALPHA_FACTOR) for (o, d, t, a) in trips]
    s2 = solve_stage2(tb, s1,
                      allow_contraflow=cfg["contra"],
                      allow_v2g=cfg["v2g"],
                      theta_t_over_p=D.THETA_T_OVER_P if theta is None else theta,
                      trip_requests=trips)
    return s1, s2


def main():
    tb = D.build()
    rows = []
    store = {}
    for cid in sorted(CASES):
        s1, s2 = run_case(tb, cid)
        store[cid] = (s1, s2)
        rows.append(dict(
            case=cid,
            Rf=100 * s1.failure_rate,
            travel=s2.travel_cost,
            loadloss=s2.load_loss_weighted,
            unserved=s2.unserved,
            Re=resilience_metric(tb, s1, s2),
        ))

    print("\n=== Stage I: UTN failure rate, Eq. (61) ===")
    print(f"  Case 1 (no crossing elimination): {rows[0]['Rf']:.2f} %"
          f"   [paper: 100.00 %]")
    print(f"  Case 2 (crossing elimination)   : {rows[1]['Rf']:.2f} %"
          f"   [paper:  38.53 %]")
    red = 100 * (rows[0]['Rf'] - rows[1]['Rf']) / rows[0]['Rf']
    print(f"  reduction                        : {red:.2f} %   [paper: 61.47 %]")

    print("\n=== Table III analogue: Stage II travel time cost ===")
    print("  case   travel cost    unserved trips")
    for r in rows:
        if r["case"] in (2, 3, 5):
            print(f"   {r['case']}     {r['travel']:11.4f}   {r['unserved']:8.2f}")
    if rows[1]["travel"] > 0:
        d = 100 * (rows[1]["travel"] - rows[2]["travel"]) / rows[1]["travel"]
        print(f"  contraflow (2 -> 3) reduces travel cost by {d:.2f} %"
              f"   [paper: 46.31 %]")

    print("\n=== Table IV analogue: PDN resilience metric, Eq. (62) ===")
    print("  case   Re          weighted load loss (MWh)")
    for r in rows:
        if r["case"] in (3, 4, 5):
            print(f"   {r['case']}     {r['Re']:.4f}      {r['loadloss']:.4f}")
    if rows[2]["Re"] > 0:
        d = 100 * (rows[3]["Re"] - rows[2]["Re"]) / rows[2]["Re"]
        print(f"  V2G (3 -> 4) improves resilience by {d:.2f} %   [paper: 16.70 %]")

    print("\n=== Case 4 detail ===")
    s1, s2 = store[4]
    print("  faulted PDN zone :", sorted(b for b, v in s1.d_pdn.items() if v))
    print("  faulted UTN nodes:", sorted(n for n, v in s1.d_utn.items() if v))
    ties = [ln for ln in tb.pdn_lines if tb.a0_pdn[ln] == 0]
    print("  tie lines closed in Stage II:", [ln for ln in ties if s2.a_pdn[ln]])
    print("  contraflow links:", sorted(s2.contraflow))
    print("  FCS net energy over Stage II (MWh, + = charging load, - = V2G):")
    dt = (D.TII_H - D.TI_H) / D.N_PERIODS
    for f in sorted(tb.fcs_bus):
        e = dt * sum(s2.pfcs[(f, t)] for t in range(1, D.N_PERIODS + 1))
        print(f"    FCS{f} (PDN bus {tb.fcs_bus[f]:2d}, UTN node {tb.fcs_node[f]:2d}): {e:+.4f}")
    s1c, s2c = store[3]
    r3 = {b for b in tb.pdn_buses
          if sum(s2c.pL[(b, t)] for t in range(1, 5)) > 1e-6}
    r4 = {b for b in tb.pdn_buses
          if sum(s2.pL[(b, t)] for t in range(1, 5)) > 1e-6}
    print("  buses shed in Case 3 but fully restored in Case 4:", sorted(r3 - r4))

    if "--sweep" in sys.argv:
        print("\n=== Table VI analogue: theta_T / theta_P sensitivity (Case 4) ===")
        print("  theta_T/theta_P   weighted load loss   travel cost")
        for th in (0.0, 1e-4, 1.0, 1e4, 1e8):
            _, s2s = run_case(tb, 4, theta=th)
            print(f"   {th:<14g}  {s2s.load_loss_weighted:16.4f}   {s2s.travel_cost:11.4f}")

    if "--plots" in sys.argv:
        from ptn.plots import make_plots
        make_plots(tb, store)
        print("\n  figures written to out/")


if __name__ == "__main__":
    main()
