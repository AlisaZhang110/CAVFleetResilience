"""
Minimal algebraic-modelling layer on top of scipy.optimize.milp (HiGHS).

Exists so that stage1.py / stage2.py can be written to look like the paper's
equations instead of like sparse-matrix assembly code.

Usage
-----
    m = Model()
    x = m.var("x", lb=0, ub=1, integer=True)
    m.constr({x: 1.0, y: -1.0}, "<=", 0.0, name="eq9_upper")
    m.minimize({x: 3.0})
    sol = m.solve()
    sol[x]
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

INF = np.inf


class Model:
    def __init__(self, name: str = "model"):
        self.name = name
        self._lb: list[float] = []
        self._ub: list[float] = []
        self._int: list[int] = []
        self._names: list[str] = []
        # constraints stored as (row_dict, lo, hi, name)
        self._rows: list[tuple[dict, float, float, str]] = []
        self._obj: dict[int, float] = {}
        self._res = None

    # ---------------- variables ----------------
    def var(self, name, lb=0.0, ub=INF, integer=False) -> int:
        idx = len(self._lb)
        self._lb.append(lb)
        self._ub.append(ub)
        self._int.append(1 if integer else 0)
        self._names.append(name)
        return idx

    def binary(self, name) -> int:
        return self.var(name, 0.0, 1.0, integer=True)

    def free(self, name) -> int:
        return self.var(name, -INF, INF)

    # ---------------- constraints ----------------
    def constr(self, expr: dict, sense: str, rhs: float, name: str = ""):
        """expr: {var_index: coeff}.  sense in {'<=', '>=', '=='}."""
        expr = {k: v for k, v in expr.items() if v != 0.0}
        if sense == "<=":
            lo, hi = -INF, float(rhs)
        elif sense == ">=":
            lo, hi = float(rhs), INF
        elif sense in ("==", "="):
            lo = hi = float(rhs)
        else:
            raise ValueError(sense)
        if not expr:
            # constant constraint: check feasibility now
            if not (lo - 1e-9 <= 0.0 <= hi + 1e-9):
                raise ValueError(f"infeasible constant constraint {name}")
            return
        self._rows.append((expr, lo, hi, name))

    def range_constr(self, expr: dict, lo: float, hi: float, name: str = ""):
        expr = {k: v for k, v in expr.items() if v != 0.0}
        if expr:
            self._rows.append((expr, float(lo), float(hi), name))

    # ---------------- objective ----------------
    def minimize(self, expr: dict):
        self._obj = dict(expr)

    def add_obj(self, expr: dict):
        for k, v in expr.items():
            self._obj[k] = self._obj.get(k, 0.0) + v

    # ---------------- solve ----------------
    def solve(self, time_limit=300.0, mip_gap=1e-4, verbose=False):
        n = len(self._lb)
        c = np.zeros(n)
        for k, v in self._obj.items():
            c[k] += v

        rows, cols, vals, lo, hi = [], [], [], [], []
        for i, (expr, l, h, _) in enumerate(self._rows):
            for k, v in expr.items():
                rows.append(i)
                cols.append(k)
                vals.append(v)
            lo.append(l)
            hi.append(h)
        A = coo_matrix((vals, (rows, cols)), shape=(len(self._rows), n)).tocsr()
        cons = LinearConstraint(A, np.array(lo), np.array(hi))

        res = milp(
            c=c,
            constraints=cons,
            bounds=Bounds(np.array(self._lb), np.array(self._ub)),
            integrality=np.array(self._int),
            options=dict(time_limit=time_limit, mip_rel_gap=mip_gap, disp=verbose),
        )
        self._res = res
        if res.x is None:
            raise RuntimeError(f"[{self.name}] no solution: {res.message}")
        return Solution(res.x, res.fun, self._names, res.message)

    @property
    def n_vars(self):
        return len(self._lb)

    @property
    def n_int(self):
        return int(sum(self._int))

    @property
    def n_constr(self):
        return len(self._rows)


class Solution:
    def __init__(self, x, obj, names, message):
        self.x = np.asarray(x)
        self.obj = obj
        self.names = names
        self.message = message

    def __getitem__(self, idx):
        return float(self.x[idx])

    def get(self, idx, default=0.0):
        return float(self.x[idx]) if idx is not None else default

    def b(self, idx, tol=0.5):
        """Read a binary variable as a bool."""
        return self.x[idx] > tol
