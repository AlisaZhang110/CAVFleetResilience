"""Figures analogous to Figs. 7, 9 and 10 of Liu et al. (2025)."""
from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

from . import data as D

OUT = "out"


def _utn_layout(tb):
    g = nx.Graph()
    g.add_nodes_from(tb.utn_nodes)
    g.add_edges_from(D.utn_undirected(tb))
    return g, nx.spring_layout(g, seed=7, iterations=400)


def fig_utn(tb, store):
    g, pos = _utn_layout(tb)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for ax, cid, title in zip(axes, (1, 2),
                              ("Case 1 - no crossing elimination",
                               "Case 2 - crossing elimination")):
        s1 = store[cid][0]
        faulted = {n for n, v in s1.d_utn.items() if v}
        alive = [(i, j) for (i, j) in g.edges
                 if s1.qstar[(i, j)] > 1e-6 or s1.qstar[(j, i)] > 1e-6]
        dead = [e for e in g.edges if e not in alive]
        nx.draw_networkx_edges(g, pos, edgelist=dead, ax=ax,
                               edge_color="#cccccc", style="dashed", width=1.2)
        nx.draw_networkx_edges(g, pos, edgelist=alive, ax=ax,
                               edge_color="#2c6fbb", width=2.4)
        cols = ["#d94040" if n in faulted else "#f2f2f2" for n in g.nodes]
        nx.draw_networkx_nodes(g, pos, ax=ax, node_color=cols,
                               edgecolors="#333333", node_size=520)
        fcs = [n for n in tb.fcs_node.values()]
        nx.draw_networkx_nodes(g, pos, nodelist=fcs, ax=ax, node_shape="s",
                               node_color="none", edgecolors="#e08a00",
                               linewidths=2.4, node_size=760)
        nx.draw_networkx_labels(g, pos, ax=ax, font_size=10)
        ax.set_title(f"{title}\nfailure rate Rf = {100*s1.failure_rate:.2f} %")
        ax.axis("off")
    fig.suptitle("UTN topology after Stage I  (square outline = FCS node, red = faulted)")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_utn_stage1.png", dpi=150)
    plt.close(fig)


def fig_v2g(tb, store):
    periods = list(range(1, D.N_PERIODS + 1))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    for ax, cid in zip(axes, (4, 5)):
        s2 = store[cid][1]
        w = 0.2
        for k, t in enumerate(periods):
            vals = [-s2.pfcs[(f, t)] for f in sorted(tb.fcs_bus)]
            ax.bar([i + (k - 1.5) * w for i in range(len(vals))], vals,
                   width=w, label=f"period {t}")
        ax.set_xticks(range(len(tb.fcs_bus)))
        ax.set_xticklabels([f"FCS{f}" for f in sorted(tb.fcs_bus)])
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(f"Case {cid}")
        ax.set_xlabel("fast charging station")
    axes[0].set_ylabel("V2G power (MW)   [+ = discharging]")
    axes[0].legend(fontsize=8)
    fig.suptitle("V2G power of AEVs")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_v2g.png", dpi=150)
    plt.close(fig)


def fig_loss(tb, store):
    periods = list(range(1, D.N_PERIODS + 1))
    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True, sharey=True)
    for ax, cid in zip(axes, (3, 4, 5)):
        s2 = store[cid][1]
        w = 0.2
        for k, t in enumerate(periods):
            vals = [s2.pL[(b, t)] for b in tb.pdn_buses]
            ax.bar([i + (k - 1.5) * w for i in range(len(vals))], vals,
                   width=w, label=f"period {t}")
        ax.set_ylabel("power loss (MW)")
        ax.set_title(f"Case {cid}", loc="left", fontsize=10)
    axes[-1].set_xticks(range(len(tb.pdn_buses)))
    axes[-1].set_xticklabels(tb.pdn_buses, fontsize=7)
    axes[-1].set_xlabel("power bus")
    axes[0].legend(fontsize=8, ncol=4)
    fig.suptitle("PDN power loss in Stage II")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_pdn_loss.png", dpi=150)
    plt.close(fig)


def make_plots(tb, store):
    os.makedirs(OUT, exist_ok=True)
    fig_utn(tb, store)
    fig_v2g(tb, store)
    fig_loss(tb, store)
