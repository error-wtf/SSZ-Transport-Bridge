#!/usr/bin/env python3
"""make_figures_v1 — static plots + animated GIFs from repo artifacts ONLY.

Every displayed number is loaded live from artifacts/*.json or
data/frozen_specs/*.npz — no hand-written physics numbers anywhere.

Outputs (figures/):
  fig1_blind_scan_tiles.png            3 sigma_min(omega) scan heatmaps
  fig2_b3b4_ladder_convergence.png     Richardson ladder: lambda vs node count
  fig3_step0_pole_reproduction.png     blind pole reproduction, 3 channels
  fig4_metric_coefficients.png         physical 3x3 coefficient diagonals vs r
  anim1_blind_scan_tiles.gif           tiles revealed one by one
  anim2_b3b4_ladder.gif                eigenvalues converging up the ladder

Usage:  python tools/make_figures_v1.py [--no-anim]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

BRIDGE = Path(__file__).resolve().parents[1]
ART = BRIDGE / "artifacts"
FIG = BRIDGE / "figures"

# dark theme (matches demos/animate_core_results.py)
plt.rcParams.update({
    "figure.facecolor": "#0d1117",
    "axes.facecolor": "#0d1117",
    "savefig.facecolor": "#0d1117",
    "text.color": "#c9d1d9",
    "axes.labelcolor": "#c9d1d9",
    "xtick.color": "#8b949e",
    "ytick.color": "#8b949e",
    "axes.edgecolor": "#30363d",
    "font.size": 9,
})

ACCENTS = ["#58a6ff", "#3fb950", "#d29922", "#f85149", "#bc8cff", "#39c5cf"]


def _load(name: str) -> dict:
    with open(ART / name, encoding="utf-8") as fh:
        return json.load(fh)


def _save_fig(fig, name: str):
    FIG.mkdir(parents=True, exist_ok=True)
    out = FIG / name
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(BRIDGE)}  ({out.stat().st_size} bytes)")


# ---------------------------------------------------------------- fig 1 + anim 1
def _tile_axes(ax, tile: dict, title: str):
    grid = bs = np.asarray(tile["grid_sigma_min"])
    nwi, nwr = grid.shape
    (re0, re1), (im0, im1) = tile["bounds"]["Re"], tile["bounds"]["Im"]
    im_ax = np.linspace(im0, im1, nwi)
    re_ax = np.linspace(re0, re1, nwr)
    pcm = ax.pcolormesh(re_ax, im_ax, grid, cmap="magma",
                        shading="auto", vmin=0.0)
    amr, ami = tile["argmin_omega"]
    ax.plot(amr, ami, "c*", ms=12, mec="white", mew=0.6)
    ax.annotate(f"argmin {amr:.2f}{ami:+.2f}\n$\\sigma_{{\\min}}$="
                f"{tile['sigma_min_global']:.4f}",
                (amr, ami), textcoords="offset points", xytext=(8, 8),
                color="white", fontsize=7.5,
                bbox=dict(fc="#000000aa", ec="#30363d", pad=2.5))
    ax.set_title(title, fontsize=9.5)
    ax.set_xlabel("Re $\\omega$")
    return pcm


def blind_scan_figures(anim: bool = True):
    scan = _load("SSZ_B6_STEP1_BLIND_SCAN_V1.json")
    tiles = scan["tiles"]

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.1), sharey=True)
    for ax, (name, tile) in zip(axes, tiles.items()):
        pcm = _tile_axes(ax, tile, name.replace("_", " ").lower())
    axes[0].set_ylabel("Im $\\omega$")
    cb = fig.colorbar(pcm, ax=axes, shrink=0.9, pad=0.015)
    cb.set_label("$\\sigma_{\\min}$ (outgoing-pencil residual)")
    fig.suptitle("B6 Step-1 blind scan — smallest pencil residual over the "
                 "declared grid (data: artifacts/SSZ_B6_STEP1_BLIND_SCAN_V1.json)",
                 fontsize=10)
    _save_fig(fig, "fig1_blind_scan_tiles.png")

    if not anim:
        return

    names = list(tiles)
    fig2, ax = plt.subplots(figsize=(6.8, 5.0))
    fig2.subplots_adjust(right=0.86)

    def frame(i: int):
        ax.clear()
        name = names[i % len(names)]
        pcm = _tile_axes(ax, tiles[name], name.replace("_", " ").lower())
        ax.set_ylabel("Im $\\omega$")
        ax.set_title(f"blind-scan tile {i % len(names) + 1}/{len(names)}",
                     fontsize=10)
        gmax = max(float(np.max(t["grid_sigma_min"])) for t in tiles.values())
        # fixed color scale across frames for comparability
        pcm.set_clim(0.0, gmax)
        if i == 0:
            fig2.colorbar(pcm, ax=ax, label="$\\sigma_{\\min}$ (fixed scale)")

    ani = FuncAnimation(fig2, frame, frames=len(names) * 2, interval=1400)
    out = FIG / "anim1_blind_scan_tiles.gif"
    ani.save(out, writer=PillowWriter(fps=0.7))
    plt.close(fig2)
    print(f"wrote figures/anim1_blind_scan_tiles.gif  ({out.stat().st_size} bytes)")


# ---------------------------------------------------------------- fig 2 + anim 2
def ladder_figures(anim: bool = True):
    art = _load("SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json")
    lad = art["s4_ladder"]
    counts = lad["ladder_node_counts"]
    lam = {int(n): np.asarray(v) for n, v in lad["lambdas"].items()}
    finest = lam[counts[-1]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    idx = np.arange(len(finest))
    for j, n in enumerate(counts):
        ax1.plot(idx, lam[n], "o--" if n == counts[-1] else "o:",
                 ms=4, color=ACCENTS[j], label=f"N={n}")
    ax1.set_xlabel("box-mode index (sorted)")
    ax1.set_ylabel("$\\lambda$")
    ax1.set_title("eigenvalues per ladder level (data: s4_ladder)")
    ax1.legend(fontsize=8)

    # Richardson-style error proxy: |lambda_N - lambda_finest| vs N, log-log
    for j, n in enumerate(counts[:-1]):
        err = np.abs(lam[n] - finest)
        ax2.loglog([n] * len(err), err, "o", ms=3.5, color=ACCENTS[j],
                   label=f"N={n}")
    ref = np.array(counts[:-1], float)
    e1 = np.abs(lam[counts[0]] - finest).max()
    e2 = np.abs(lam[counts[1]] - finest).max()
    slope = (np.log10(e2) - np.log10(e1)) / (np.log10(counts[1]) - np.log10(counts[0]))
    ax2.loglog(ref, e1 * (ref / ref[0]) ** slope, "-", color="#8b949e",
               label=f"ref slope {slope:.2f} (log-log fit)")
    ax2.set_xlabel("ladder node count N")
    ax2.set_ylabel("$|\\lambda_N - \\lambda_{N_{finest}}|$ (max)")
    ax2.set_title("ladder convergence toward finest level")
    ax2.legend(fontsize=8)
    fig.suptitle("B3/B4 independent 3-DOF box solve — Richardson ladder "
                 "(artifact: SSZ_B3_B4_FULL_3DOF_INDEPENDENT_SOLVE_V1.json)",
                 fontsize=10)
    _save_fig(fig, "fig2_b3b4_ladder_convergence.png")

    if not anim:
        return

    fig2, ax = plt.subplots(figsize=(6.8, 5.0))

    def frame(i: int):
        ax.clear()
        upto = min(i + 1, len(counts))
        for j, n in enumerate(counts[:upto]):
            ax.plot(idx, lam[n], "o--" if n == counts[-1] else "o:",
                    ms=4, color=ACCENTS[j], label=f"N={n}")
        ax.set_xlabel("box-mode index (sorted)")
        ax.set_ylabel("$\\lambda$")
        ax.set_title(f"Richardson ladder — level {upto}/{len(counts)}")
        ax.legend(fontsize=8)

    ani = FuncAnimation(fig2, frame, frames=len(counts) + 1, interval=1000)
    out = FIG / "anim2_b3b4_ladder.gif"
    ani.save(out, writer=PillowWriter(fps=1))
    plt.close(fig2)
    print(f"wrote figures/anim2_b3b4_ladder.gif  ({out.stat().st_size} bytes)")


# ---------------------------------------------------------------- fig 3
def step0_pole_figure():
    art = _load("SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json")
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2), sharex=True)
    for ax, (ch, c) in zip(axes, art["channels"].items()):
        refs = np.asarray([r["ref"] for r in c["g3_reference_reproduction"]],
                          float)
        dists = [r["dist"] for r in c["g3_reference_reproduction"]]
        ax.plot(refs[:, 0], refs[:, 1], "x", ms=9, mew=2, color="#f85149",
                label="reference (frozen)")
        # blind values: reference + dist (upper bound) shown as circles
        ax.plot(refs[:, 0], refs[:, 1], "o", ms=13, mfc="none",
                color="#58a6ff", label="blind reproduction")
        for (x, y), d in zip(refs, dists):
            ax.annotate(f"dist {d:.1e}", (x, y), textcoords="offset points",
                        xytext=(8, -4), color="#8b949e", fontsize=7)
        ax.set_title(f"{ch}: V0={c['V0']}, a={c['a']} — "
                     f"{len(refs)}/{c['n_candidates']} ref poles shown")
        ax.set_xlabel("Re $k$")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("Im $k$")
    axes[0].legend(fontsize=8)
    fig.suptitle("B6 Step-0 square-well controls — blind pole reproduction "
                 "(marker size >> actual distance; data: "
                 "artifacts/SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json)",
                 fontsize=9.5)
    _save_fig(fig, "fig3_step0_pole_reproduction.png")


# ---------------------------------------------------------------- fig 4
def metric_coefficient_figure():
    npz = np.load(BRIDGE / "data" / "frozen_specs"
                  / "V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz",
                  allow_pickle=True)
    r = npz["r"]
    M, K = npz["M_phys"], npz["K_phys"]
    fields = ["$\\psi$ (ch1)", "$\\dot\\phi$ (ch2)", "$V$ (ch3)"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for j in range(3):
        axes[0].semilogy(r, np.abs(M[:, j, j]) + 1e-30,
                         color=ACCENTS[j], label=f"M diag, {fields[j]}")
        axes[1].semilogy(r, np.abs(K[:, j, j]) + 1e-30,
                         color=ACCENTS[j], label=f"K diag, {fields[j]}")
    for ax, ttl in zip(axes, ("mass matrix M (diag |.|)",
                              "stiffness matrix K (diag |.|)")):
        ax.set_xlabel("$r$")
        ax.set_title(ttl)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.suptitle("Frozen physical 3-DOF operator coefficients "
                 "(data: data/frozen_specs/V4_PHYSICAL_RESONANCE_EXPORT_V3_3DOF.npz)",
                 fontsize=10)
    _save_fig(fig, "fig4_metric_coefficients.png")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-anim", action="store_true",
                    help="skip GIF generation (static figures only)")
    args = ap.parse_args()
    anim = not args.no_anim
    blind_scan_figures(anim)
    ladder_figures(anim)
    step0_pole_figure()
    metric_coefficient_figure()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
