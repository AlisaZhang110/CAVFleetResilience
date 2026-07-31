# Replication: Coordinated Topology Reconfiguration for Resilience Enhancement in Coupled Power–Transportation Networks Based on Autonomous EV Fleets

Liu, Xue, Zhou, Han, Chang, Li, Su, Shahidehpour & Sun, *IEEE Trans. Transportation Electrification* **11**(4):10062–10075, Aug. 2025. DOI 10.1109/TTE.2025.3560376

```
python3 run_cases.py                 # Cases 1-5, Tables III/IV analogues
python3 run_cases.py --sweep         # + theta_T/theta_P sensitivity (Table VI)
python3 run_cases.py --plots         # + figures in out/
```

Solver: HiGHS via `scipy.optimize.milp`. The paper uses GUROBI through MATLAB; no
commercial solver is needed at this scale (Stage II ≈ 3.5k variables, 150 binaries,
solves in well under a second).

```
ptn/lp.py        thin algebraic layer over scipy.optimize.milp
ptn/data.py      P33T12 testbed, every value tagged [PAPER] / [STANDARD] / [ASSUMED]
ptn/stage1.py    Eqs. (1)-(25)   fault isolation
ptn/stage2.py    Eqs. (26)-(60)  service restoration + AEV fleet + V2G
ptn/plots.py     Figs. 7 / 9 / 10 analogues
run_cases.py     Table II case matrix and the reported metrics
```

---

## 1. What the paper does

A two-stage centralised recovery model for a power distribution network (PDN)
coupled to an urban transportation network (UTN) through fast charging stations
(FCSs), after an extreme-weather fault.

**Stage I — fault isolation (15 min, single period).** Faulted zones grow along
still-closed lines/links (Eqs. 1–4). A faulted PDN bus hosting an FCS also faults
the corresponding UTN node (Eq. 5) — this is the only cross-network propagation
direction modelled. Two actions run simultaneously: PDN *line switching* (open
sectionalising switches to bound the zone) and UTN *crossing elimination* (place
roadblocks so congestion around dead FCSs does not spread). Objective (25):
minimise weighted load curtailment, maximise passable road capacity. Stage I can
only **open** lines, never close them (Eq. 9) — restoration is deferred.

**Stage II — service restoration (4 h, 4 periods).** Takes Stage I's topology as
input. PDN line switching may now **close** tie lines (Eq. 26) subject to radiality
(27–28); UTN *link contraflow* reverses lanes to double capacity in the needed
direction (44–46). Simultaneously an autonomous EV fleet is dispatched on a
time–space–state-of-charge expanded graph (47–48): occupied AEVs serve exogenous
trip requests, unoccupied AEVs rebalance and either charge or discharge at FCSs.
Discharging is the V2G channel that back-feeds the PDN (Eq. 57), which is what
lets an FCS act as a generator inside an island. Objective (60): weighted sum of
PDN load loss and system-optimal (Wardrop-second-principle) AEV travel cost.

**Claimed results.** UTN failure rate 100% → 38.53% from crossing elimination;
travel cost −46.31% from contraflow; PDN resilience metric 0.7252 → 0.8463
(+16.70%) from V2G; θ_T/θ_P = 1 minimises both objective terms.

---

## 2. What is and is not replicable

**Not recoverable from the document.** The PDN load vector (Fig. 6 only, and its
peaks near 0.85 MW do *not* match the standard 33-bus feeder's 0.42 MW, so the
loads were scaled or replaced); the UTN adjacency (Fig. 5b only); link capacities
Q⁰ and free-flow times t⁰; β; the load weights ω_i and the first/second-level bus
sets; the DG buses and ratings; which lines carry sectionalising switches; the
fleet size, initial distribution and SoC discretisation; the trip request set
beyond the single 1→12 example; and the fault set beyond the one named line.
**Exact numerical replication of Tables III–VI is therefore impossible**, and any
claim to have matched them would be an artefact of tuning.

**What this code does replicate.** The full model structure, and — on a documented
synthetic testbed — every qualitative result. Two of the reconstructions are
genuinely load-bearing checks:

* The Stage I fault-propagation model, given a fault on line 14–15 and switches on
  the neighbouring lines, produces exactly the paper's isolated zone behaviour, and
  the Case 1 UTN failure rate comes out at **exactly 100%**.
* Stage II independently selects contraflow on links **2→1, 6→2, 7→6, 12→7** —
  precisely the four the paper names in Sec. IV-B-1 — because those are the reverse
  directions along the stated evacuation route 1→2→6→7→12.

The fault scenario in `data.py` was reverse-engineered so that the paper's Sec. IV-B
narrative is internally consistent: faults on (9,10), (14,15) and (28,29) produce
faulted FCSs on the PDN buses the paper names, an island behind the first fault
(so a bus can be "restored in Case 4"), and the 30–33 island that FCS6 back-feeds
while line 18–33 opens. Only (14,15) is stated in the paper.

### Results on this testbed

| quantity | this replication | paper |
|---|---|---|
| UTN failure rate, Case 1 | 100.00 % | 100.00 % |
| UTN failure rate, Case 2 | 50.00 % | 38.53 % |
| reduction from crossing elimination | 50.00 % | 61.47 % |
| travel cost reduction from contraflow (2→3) | 49.4 % | 46.31 % |
| resilience gain from V2G (3→4) | 12.0 % | 16.70 % |
| contraflow links chosen | 2→1, 6→2, 7→6, 12→7 (+3 for a second trip) | 2→1, 6→2, 7→6, 12→7 |
| FCSs discharging in Case 4 | FCS3, FCS6 (and FCS1, FCS5) | FCS3, FCS6 |
| FCSs unable to exchange | FCS2, FCS4 | FCS2, FCS4 |
| θ_T/θ_P sweep shape | load loss ↓ with θ_P, travel cost ↓ with θ_T | same |

**One claim does not reproduce.** The paper reports Case 5 (robust) resilience
*worse* than Case 3, attributing it to extra charging demand from the larger trip
requests. Here Case 5 has identical load loss to Case 4, only higher travel cost.
The reason is structural: in the paper's Eqs. (54)/(57) only *unoccupied* AEVs
appear at FCSs, so inflating occupied demand cannot by itself raise charging load.
Getting the paper's effect requires occupied vehicles to recharge, which is not in
the printed formulation.

---

## 3. Errors and ambiguities in the source paper

Flagged because they change the implementation, not for their own sake.

| # | Location | Issue |
|---|---|---|
| DEV-1 | Eqs. (3)(4) | Propagation is written against `a_ij,0`, the *original* line status. With that form no switching action could ever bound a zone and Cases 1 and 2 would be identical. Stage II's Eqs. (29)(30) use the post-action status, so `a_ij,I` is clearly intended. Implemented that way. |
| DEV-2 | Eq. (5) | Written as an equality `d_q = d_r`, which forces UTN→PDN propagation. The surrounding text explicitly says that direction is excluded — and then says so *again* two paragraphs later while leaving the equality in place. Implemented as the one-way implication `d_r ≥ d_q`. |
| DEV-3 | Eq. (21) | `a_ij,I = e_i(j)` inverts the meaning: `a = 1` is "unobstructed", `e = 1` is "blocked". Sign corrected to `a = 1 − e`. |
| DEV-4 | Eq. (22) | `Σ e_i(k) ≤ N_c`. `N_c` is undefined in the nomenclature (`C` is the SoC set; `N_i` is the intersection degree), and as written the constraint is vacuous — you can never block more links than exist. Implemented as a budget of `N_i − 1` at healthy crossings, unrestricted at faulted ones (which is the action the text describes). |
| DEV-5 | Eq. (25) | The two terms are weighted MWh of lost load and a raw vehicle-capacity count, summed without weights. Dimensionally incoherent, and the relative scaling silently determines the answer. Stage II adds θ_P/θ_T for exactly this reason; Stage I has no equivalent. Implemented with the UTN term normalised to a failure fraction and a documented scale factor. |
| DEV-6 | Eq. (47) | `L_R` admits only movement links, so a vehicle must traverse a link every period and cannot dwell. Infeasible for any horizon longer than the shortest path. Zero-cost waiting links added. |
| DEV-7 | Eq. (49) | Two problems. (a) The BPR function is truncated to `t* = β t⁰ (f/Q)` — the standard `t⁰(1 + 0.15(v/c)⁴)` has lost both its constant term and its exponent, so free-flow time is zero at zero flow. (b) β is simultaneously the "economic value of time" per the nomenclature, so (49) returns money while being used as a time. Consequently `f · t*(f)` in Eq. (58) is **quadratic**, meaning the printed Stage II problem is an MIQP, not the MILP implied by "convex mixed-integer model … handled centrally by the solver". Kept as printed; the convex quadratic is handled by tangent-plane outer linearisation. |
| DEV-8 | Eqs. (54)(49) | Only unoccupied flow `f_e` is charged against link capacity and travel cost. Occupied AEVs are invisible to congestion — yet the entire contraflow result is about accelerating an occupied trip. Occupied flows counted here too. |
| DEV-9 | Eq. (62) | Uses `B` as a bus set; `B` is defined in the nomenclature as the set of blocked *links*. Also, the text calls it "the ratio of total restored load expectations", but the formula is `1 − loss/total`. |
| — | Sec. V | The conclusion states Stage I reduces "PDN load losses by 61.47%". 61.47% is the reduction in the **UTN failure rate** (Case 1 → Case 2, Sec. IV-B-1). No PDN load-loss figure of that size appears anywhere in the paper. |
| — | Tables V, VI | Identical captions ("Comparison of Computational Performance in Cases 4 and 6"), though Table VI is the weight sensitivity and only five cases are ever defined. |
| — | Eq. (46) | `Σ c_ij ≤ N_l` where `N_l` is the total number of links — vacuous, as with Eq. (22). Presumably intended as a reversal budget. |
| — | Appendix | With a box uncertainty set and no budget-of-uncertainty parameter, the dual reformulation (66)–(69) collapses to substituting `ᾱ_o`. The RO machinery adds nothing beyond evaluating the model at the worst-case demand; Case 5 is implemented that way. |
| — | Eqs. (38)-(41) | Reference `d_i,s,II`, but no Stage II constraint ever defines it — Eqs. (29)(30) are written in `d_i,s,I`. Treated as the Stage I zone carried forward. |

Separately worth checking against the sources rather than taking on trust: the
claim that cutting off a faulted FCS "will not affect the operation of the PDN" is
cited to [30] (Liao, Taiebat & Xu, *Appl. Energy* 2021), which is an economic and
environmental viability study of shared AEV fleets — not an obvious support for a
protection-coordination claim.

---

## 4. Modelling choices we may want to change

* `data.py :: PDN_FAULTS`, `UTN_FAULTS` — the scenario. Only (14,15) is the paper's.
* `data.py :: SWITCHED_LINES` — controls how tightly zones can be bounded; this is
  the single most influential assumed input for the Stage I numbers.
* `data.py :: LOAD_SCALE` — set to 2.0 to match Fig. 6's magnitudes rather than the
  textbook 33-bus loads. The V2G resilience gain is sensitive to this.
* `data.py :: LAMBDA_CAP_STAGE1` — the DEV-5 scale factor.
* `data.py :: SOC_LEVELS`, `FLEET_INIT`, `TRIP_REQUESTS` — fleet model granularity.
* `stage2.py :: n_tangents` — accuracy of the DEV-7 linearisation (outer
  approximation, so it slightly *under*-estimates travel cost; raise for tightness).
