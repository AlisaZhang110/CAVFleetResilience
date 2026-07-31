"""
P33T12 testbed: 33-bus PDN coupled to a 12-node UTN.

PROVENANCE TAGGING
------------------
Every quantity below carries one of three tags in its comment:

  [PAPER]    Stated numerically in Liu et al. (2025), IEEE TTE 11(4):10062-10075.
  [STANDARD] Taken from the standard IEEE 33-bus test feeder (Baran & Wu, 1989),
             which the paper cites as its PDN but does not tabulate.
  [ASSUMED]  NOT recoverable from the paper. Chosen here to be self-consistent
             with the qualitative narrative in Section IV-B. Change freely.

The paper does not tabulate its PDN loads (only Fig. 6), its UTN topology
(only Fig. 5b), its link capacities, free-flow times, trip requests, fleet
size, load weights, or fault set (only Fig. 5 pictograms plus the single
sentence naming Line 14-15). Exact numerical replication is therefore not
possible from the document alone. See README.md, "What is and is not
replicable".
"""
from __future__ import annotations

from dataclasses import dataclass, field

# ============================================================================
# 1. PDN: IEEE 33-bus
# ============================================================================

SBASE = 10.0      # MVA          [STANDARD]
VBASE = 12.66     # kV           [STANDARD]
ZBASE = VBASE ** 2 / SBASE  # 16.0276 ohm

# (from, to, R[ohm], X[ohm])                                        [STANDARD]
PDN_BRANCHES = [
    (1, 2, 0.0922, 0.0477), (2, 3, 0.4930, 0.2511), (3, 4, 0.3660, 0.1864),
    (4, 5, 0.3811, 0.1941), (5, 6, 0.8190, 0.7070), (6, 7, 0.1872, 0.6188),
    (7, 8, 0.7114, 0.2351), (8, 9, 1.0300, 0.7400), (9, 10, 1.0440, 0.7400),
    (10, 11, 0.1966, 0.0650), (11, 12, 0.3744, 0.1238), (12, 13, 1.4680, 1.1550),
    (13, 14, 0.5416, 0.7129), (14, 15, 0.5910, 0.5260), (15, 16, 0.7463, 0.5450),
    (16, 17, 1.2890, 1.7210), (17, 18, 0.7320, 0.5740), (2, 19, 0.1640, 0.1565),
    (19, 20, 1.5042, 1.3554), (20, 21, 0.4095, 0.4784), (21, 22, 0.7089, 0.9373),
    (3, 23, 0.4512, 0.3083), (23, 24, 0.8980, 0.7091), (24, 25, 0.8960, 0.7011),
    (6, 26, 0.2030, 0.1034), (26, 27, 0.2842, 0.1447), (27, 28, 1.0590, 0.9337),
    (28, 29, 0.8042, 0.7006), (29, 30, 0.5075, 0.2585), (30, 31, 0.9744, 0.9630),
    (31, 32, 0.3105, 0.3619), (32, 33, 0.3410, 0.5302),
]

# Normally-open tie lines (from, to, R, X)                          [STANDARD]
PDN_TIES = [
    (8, 21, 2.0, 2.0), (9, 15, 2.0, 2.0), (12, 22, 2.0, 2.0),
    (18, 33, 0.5, 0.5), (25, 29, 0.5, 0.5),
]

# bus -> (P[MW], Q[MVar])                                           [STANDARD]
# NOTE: the paper's Fig. 6 shows peaks near 0.85 MW at several buses, which
# does NOT match this standard set (peak 0.42 MW). The paper evidently scaled
# or replaced the loads.  Use LOAD_SCALE to explore.
PDN_LOADS = {
    2: (0.100, 0.060), 3: (0.090, 0.040), 4: (0.120, 0.080), 5: (0.060, 0.030),
    6: (0.060, 0.020), 7: (0.200, 0.100), 8: (0.200, 0.100), 9: (0.060, 0.020),
    10: (0.060, 0.020), 11: (0.045, 0.030), 12: (0.060, 0.035), 13: (0.060, 0.035),
    14: (0.120, 0.080), 15: (0.060, 0.010), 16: (0.060, 0.020), 17: (0.060, 0.020),
    18: (0.090, 0.040), 19: (0.090, 0.040), 20: (0.090, 0.040), 21: (0.090, 0.040),
    22: (0.090, 0.040), 23: (0.090, 0.050), 24: (0.420, 0.200), 25: (0.420, 0.200),
    26: (0.060, 0.025), 27: (0.060, 0.025), 28: (0.060, 0.020), 29: (0.120, 0.070),
    30: (0.200, 0.600), 31: (0.150, 0.070), 32: (0.210, 0.100), 33: (0.060, 0.040),
}
LOAD_SCALE = 2.0  # [ASSUMED] Fig. 6 of the paper peaks near 0.85 MW, roughly 2x the
                  #           standard 33-bus peak of 0.42 MW, so the paper evidently
                  #           scaled or replaced the loads.

SUBSTATION_BUS = 1                                       # [STANDARD]
# Distributed generators.  Fig. 5(a) shows generator pictograms but the buses
# are not legible.  Bus 18 is required for the Case 3 narrative (tie 18-33
# closed to feed the 30-33 island); bus 22 backs the tie 12-22 restoration.
DG_BUSES = {18: (0.60, 0.40), 22: (0.60, 0.40)}          # bus -> (Pmax, Qmax) MW  [ASSUMED]
SUBSTATION_CAP = (10.0, 8.0)                             # [ASSUMED]

# Load priority weights w_i. Fig. 5(a) marks "First level load" and
# "Second level load" but the buses are not legible.
FIRST_LEVEL = {7, 8, 24, 25, 30, 32}                     # [ASSUMED]
SECOND_LEVEL = {11, 14, 18, 22, 33}                      # [ASSUMED]
W_FIRST, W_SECOND, W_NORMAL = 10.0, 5.0, 1.0             # [ASSUMED]

VMIN, VMAX, VREF = 0.95, 1.05, 1.00                      # [ASSUMED, conventional]
SMAX_LINE = 4.0                                          # MVA  [ASSUMED]
SMAX_TIE = 2.0                                           # MVA  [ASSUMED]

# Lines carrying a remote-controlled sectionalizing switch (k_ij,0 = 1).
# All ties are switchable by definition.  Line 14-15 deliberately has NO
# switch, otherwise the fault it carries would be isolated with zero
# consequence and the paper's Case 1 narrative could not arise.
SWITCHED_LINES = {                                       # [ASSUMED]
    (2, 3), (6, 7), (7, 8), (8, 9), (10, 11), (13, 14), (15, 16),
    (19, 20), (21, 22), (23, 24), (26, 27), (27, 28), (28, 29),
    (29, 30), (30, 31), (32, 33),
}

# ============================================================================
# 2. UTN: 12 nodes
# ============================================================================

# Undirected edges -> 40 directed links, matching the ~40 link labels visible
# in Fig. 5(b).  The actual adjacency in Fig. 5(b) is not legible; this graph
# is constructed to contain the one route the paper names explicitly,
# 1 -> 2 -> 6 -> 7 -> 12, and to give the FCS nodes plausible degree.
UTN_EDGES = [                                            # [ASSUMED]
    (1, 2), (1, 3), (1, 9), (2, 5), (2, 6), (3, 4), (3, 10), (4, 5),
    (4, 9), (4, 11), (5, 6), (5, 8), (6, 7), (7, 8), (7, 12), (8, 9),
    (9, 10), (9, 12), (10, 11), (11, 12),
]

Q0_LINK = 10.0       # veh / period / direction   [ASSUMED]
T0_LINK = 1.0        # free-flow travel time, normalised to one period [ASSUMED]
BETA_VOT = 20.0      # beta, "economic value of time" ($/h)  [ASSUMED]
ZETA_LINK = 1.0      # zeta_ij, UTN link weight in Eq. (61)  [ASSUMED]

# ============================================================================
# 3. PDN <-> UTN coupling: the six fast charging stations
# ============================================================================
# The paper states (Sec. IV-B-1) that PDN buses 10 and 15 are de-energised and
# that this faults two FCS nodes in the UTN; and (Sec. IV-B-2) that "FCSs 2 and
# 4 are unable to exchange electricity between AEVs and PDN".  So
#   FCS2 <-> PDN bus 15   and   FCS4 <-> PDN bus 10   (either way round).
# Their UTN nodes must NOT lie on the route 1->2->6->7->12, because the paper
# says AEVs serve the 1->12 request along that route "while avoiding the fault
# area".  Nodes 9 and 10 satisfy that.  The remaining four are placed so that
# the rest of the narrative holds:
#   FCS6 must sit in the 30-33 island (it becomes a generation bus there);
#   FCS3 must sit somewhere that needs to discharge (bus 12, in the 11-13 island);
#   FCS1 and FCS5 must sit on healthy feeder sections (the paper has them charging).
FCS = {                                                  # name -> (pdn_bus, utn_node)
    1: (4, 3),    # [ASSUMED]
    2: (15, 10),  # [PAPER-CONSTRAINED on the PDN side]
    3: (12, 5),   # [ASSUMED]
    4: (10, 9),   # [PAPER-CONSTRAINED on the PDN side]
    5: (25, 8),   # [ASSUMED]
    6: (33, 6),   # [PAPER-CONSTRAINED to the 30-33 island]
}
FCS_CAP = 8.0        # max simultaneous vehicles at one FCS  [ASSUMED]

# ============================================================================
# 4. Fault scenario
# ============================================================================
# Only Line (14,15) is named in the paper.  Fig. 5(a) shows roughly six fault
# pictograms whose positions are not legible.  The two extra PDN faults below
# are chosen because together they reproduce every qualitative claim in
# Sec. IV-B:
#   (9,10)  -> faults bus 10 => FCS4 down; islands buses 11-13 behind it,
#              which is how bus 11 comes to be "restored in Case 4"
#   (14,15) -> faults bus 15 => FCS2 down                        [PAPER]
#   (28,29) -> islands buses 30-33, which is the "curtailed area, i.e., from
#              Buses 30 to 33" that FCS6 back-feeds in Case 4
PDN_FAULTS = [(9, 10), (14, 15), (28, 29)]               # [(14,15) PAPER; rest ASSUMED]
UTN_FAULTS = [(4, 9)]                                    # the "Damaged road" of Fig. 5(b)  [ASSUMED]

# ============================================================================
# 5. Time, fleet, demand
# ============================================================================
T0_H, TI_H, TII_H = 0.20, 0.27, 4.27                     # hours  [PAPER]
N_PERIODS = 4                                            # [PAPER] Stage II = 4 h / 4 periods
N_TIME_NODES = N_PERIODS + 1                             # 5 time layers, 4 transitions

SOC_LEVELS = 6                                           # [ASSUMED]
SOC_STEP_MW = 0.100                                      # 100 kW per unit AEV  [PAPER]
TRAVEL_SOC_COST = 1                                      # SoC levels per link traversal  [ASSUMED]

# Initial fleet: (utn_node, soc, count)                   [ASSUMED]
FLEET_INIT = [
    (1, 6, 14), (3, 6, 10), (5, 5, 6), (6, 6, 6), (8, 6, 6), (12, 6, 4),
]

# Customer trip requests (origin, destination, departure time index, alpha)
# The paper names only the aggregate request 1 -> 12.       [PAPER for #1]
TRIP_REQUESTS = [
    (1, 12, 1, 10.0),   # [PAPER: route 1->2->6->7->12, "Traffic demand O->D: 10"]
    (3, 8, 1, 6.0),     # [ASSUMED]
]

THETA_T_OVER_P = 1.0                                     # [PAPER]
# Case 5 is the robust counterpart of Case 4.  With a box uncertainty set and
# no budget-of-uncertainty constraint, the dual reformulation in the paper's
# Appendix (Eqs. 63-69) collapses to simply substituting the upper bound of
# alpha_o -- so Case 5 == Case 4 evaluated at alpha_bar.
ROBUST_ALPHA_FACTOR = 1.4                                # [ASSUMED]
LAMBDA_CAP_STAGE1 = 1.0   # scaling on the UTN term of Eq. (25)  [ASSUMED - see README]


# ============================================================================
# Derived helpers
# ============================================================================
@dataclass
class Testbed:
    pdn_buses: list = field(default_factory=list)
    pdn_lines: list = field(default_factory=list)       # (i,j)
    r: dict = field(default_factory=dict)               # pu
    x: dict = field(default_factory=dict)               # pu
    smax: dict = field(default_factory=dict)            # MVA
    a0_pdn: dict = field(default_factory=dict)          # initially closed?
    k0_pdn: dict = field(default_factory=dict)          # has switch?
    pl: dict = field(default_factory=dict)
    ql: dict = field(default_factory=dict)
    w: dict = field(default_factory=dict)
    gen: dict = field(default_factory=dict)             # bus -> (pmax,qmax)
    utn_nodes: list = field(default_factory=list)
    utn_links: list = field(default_factory=list)       # directed (i,j)
    q0: dict = field(default_factory=dict)
    t0: dict = field(default_factory=dict)
    fcs_bus: dict = field(default_factory=dict)         # fcs id -> pdn bus
    fcs_node: dict = field(default_factory=dict)        # fcs id -> utn node
    mu_pdn: dict = field(default_factory=dict)
    mu_utn: dict = field(default_factory=dict)


def build() -> Testbed:
    tb = Testbed()
    tb.pdn_buses = list(range(1, 34))

    for (i, j, r, x) in PDN_BRANCHES:
        tb.pdn_lines.append((i, j))
        tb.r[(i, j)] = r / ZBASE
        tb.x[(i, j)] = x / ZBASE
        tb.smax[(i, j)] = SMAX_LINE
        tb.a0_pdn[(i, j)] = 1
        tb.k0_pdn[(i, j)] = 1 if (i, j) in SWITCHED_LINES else 0
    for (i, j, r, x) in PDN_TIES:
        tb.pdn_lines.append((i, j))
        tb.r[(i, j)] = r / ZBASE
        tb.x[(i, j)] = x / ZBASE
        tb.smax[(i, j)] = SMAX_TIE
        tb.a0_pdn[(i, j)] = 0          # normally open
        tb.k0_pdn[(i, j)] = 1

    for b in tb.pdn_buses:
        p, q = PDN_LOADS.get(b, (0.0, 0.0))
        tb.pl[b] = p * LOAD_SCALE
        tb.ql[b] = q * LOAD_SCALE
        tb.w[b] = W_FIRST if b in FIRST_LEVEL else (
            W_SECOND if b in SECOND_LEVEL else W_NORMAL)

    tb.gen[SUBSTATION_BUS] = SUBSTATION_CAP
    for b, cap in DG_BUSES.items():
        tb.gen[b] = cap

    tb.utn_nodes = list(range(1, 13))
    for (i, j) in UTN_EDGES:
        for (a, b) in ((i, j), (j, i)):
            tb.utn_links.append((a, b))
            tb.q0[(a, b)] = Q0_LINK
            tb.t0[(a, b)] = T0_LINK

    for f, (bus, node) in FCS.items():
        tb.fcs_bus[f] = bus
        tb.fcs_node[f] = node

    for ln in tb.pdn_lines:
        tb.mu_pdn[ln] = 0
    for (i, j) in PDN_FAULTS:
        key = (i, j) if (i, j) in tb.mu_pdn else (j, i)
        tb.mu_pdn[key] = 1

    for ln in tb.utn_links:
        tb.mu_utn[ln] = 0
    for (i, j) in UTN_FAULTS:
        tb.mu_utn[(i, j)] = 1
        tb.mu_utn[(j, i)] = 1

    return tb


def utn_undirected(tb):
    seen = set()
    out = []
    for (i, j) in tb.utn_links:
        if (j, i) in seen:
            continue
        seen.add((i, j))
        out.append((i, j))
    return out


def pdn_neighbors(tb):
    """bus -> list of (line, orientation) where orientation=+1 if bus is 'from'."""
    nb = {b: [] for b in tb.pdn_buses}
    for (i, j) in tb.pdn_lines:
        nb[i].append(((i, j), +1))
        nb[j].append(((i, j), -1))
    return nb
