"""Animated GIF demos for the SSZ-Transport-Bridge core results.

Generates three dark-theme GIFs into demos/output/:

  1. b6_step0_square_well_pole_reveal.gif
     B6 step-0 analytic known-answer controls: the three frozen
     square-well channels and the six frozen reference poles revealed
     one after another, each annotated with its certified blind
     reproduction distance.
  2. bridge_contract_chain.gif
     The five structural bridge contracts lighting up in chain order;
     footer verdict loaded from the bridge certificate artifact.
  3. b6_step1_ch2_pilot_pole.gif
     B6 step-1 ch2 pilot pole: pole_scan_sign mirrored into the QNM
     sign convention, then the certification-gate captions.

Every displayed number is loaded live from artifacts/ or data/frozen_specs/
— no hand-typed values. Run:

    cd <repo> && python demos/animate_core_results.py
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyBboxPatch

BRIDGE = Path(__file__).resolve().parents[1]
ART = BRIDGE / "artifacts"
SPECS = BRIDGE / "data" / "frozen_specs"
OUT = BRIDGE / "demos" / "output"

BG = "#0e1117"
FG = "#c9d1d9"
GRID = "#30363d"
CH_COLORS = {"ch1": "#58a6ff", "ch2": "#f778ba", "ch3": "#7ee787"}
ACCENT = "#e3b341"

FPS = 2


def _load(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _dark_theme() -> None:
    plt.rcParams.update({
        "figure.facecolor": BG,
        "axes.facecolor": BG,
        "savefig.facecolor": BG,
        "text.color": FG,
        "axes.labelcolor": FG,
        "xtick.color": FG,
        "ytick.color": FG,
        "axes.edgecolor": GRID,
        "font.size": 9,
    })


def _save(anim: FuncAnimation, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    anim.save(path, writer=PillowWriter(fps=FPS), dpi=110)
    return path


def _fmt(v: float) -> str:
    return f"{v:.1e}"


# ---------------------------------------------------------------------------
# Demo 1 — B6 step 0: square-well blind pole reveal
# ---------------------------------------------------------------------------

def demo_step0_pole_reveal() -> Path:
    controls = _load(ART / "SSZ_B6_STEP0_ANALYTIC_CONTROLS_V1.json")

    gates = controls["gates"]
    n_gates = len(gates)
    n_pass = sum(1 for ok in gates.values() if ok)
    resids = [controls["channels"][ch]["g2_analytic_residual_max"]
              for ch in ("ch1", "ch2", "ch3")]
    dists = [
        entry["dist"]
        for ch in ("ch1", "ch2", "ch3")
        for entry in controls["channels"][ch]["g3_reference_reproduction"]
    ]
    neg = controls["negative_controls"]

    fig, (ax_v, ax_w) = plt.subplots(1, 2, figsize=(9.8, 4.9))
    fig.subplots_adjust(bottom=0.20)
    fig.suptitle(
        f"B6 Step 0 — analytic square-well controls: blind pole reveal "
        f"({n_pass}/{n_gates} gates PASS)",
        color=FG, fontsize=11,
    )

    # left panel: the three frozen half-line square wells
    x_max = 3.0
    for ch in ("ch1", "ch2", "ch3"):
        c = controls["channels"][ch]
        v0, a = c["V0"], c["a"]
        xsc = [0, a * 0.999, a, a, a * 1.001, x_max]
        ys = [v0, v0, v0, 0.0, 0.0, 0.0]
        ax_v.plot(xsc, ys, color=CH_COLORS[ch], lw=1.8,
                  label=f"{ch}: V0={v0:g}, a={a:g} (frozen spec)")
    ax_v.axhline(0.0, color=GRID, lw=0.8)
    ax_v.set_xlim(0, x_max)
    ax_v.set_xlabel("x")
    ax_v.set_ylabel("V(x)")
    ax_v.set_title("frozen half-line square wells (Neumann center)", fontsize=9)
    ax_v.legend(fontsize=7, loc="lower right", framealpha=0.2)
    ax_v.grid(color=GRID, lw=0.4, alpha=0.5)

    # right panel: complex plane with reference poles appearing one by one
    pole_artists = []  # (marker, label-artist)
    reveal_order = [
        (ch, i, entry)
        for ch in ("ch1", "ch2", "ch3")
        for i, entry in enumerate(controls["channels"][ch]["g3_reference_reproduction"])
    ]
    for ch, i, entry in reveal_order:
        re_w, im_w = entry["ref"]
        line, = ax_w.plot([re_w], [im_w], "o", ms=9, mfc=CH_COLORS[ch], mec="none",
                          alpha=0.0)
        txt = ax_w.annotate(
            f"{ch} ref#{i + 1}\ndist {_fmt(entry['dist'])}",
            (re_w, im_w), textcoords="offset points", xytext=(8, 6),
            color=CH_COLORS[ch], fontsize=7, alpha=0.0,
        )
        pole_artists.append((line, txt))
    ax_w.axhline(0.0, color=GRID, lw=0.8)
    ax_w.set_xlabel("Re w")
    ax_w.set_ylabel("Im w")
    ax_w.set_title("reference poles reproduced blind", fontsize=9)
    ax_w.grid(color=GRID, lw=0.4, alpha=0.5)

    footer_l1 = (
        f"G1 blind-freeze | G2 analytic residual {_fmt(resids[0])} / {_fmt(resids[1])} / "
        f"{_fmt(resids[2])} (ch1/ch2/ch3)"
    )
    footer_l2 = (
        f"G3 all {len(dists)} refs blind, dist {_fmt(min(dists))}..{_fmt(max(dists))} | "
        f"G4 neg. control: D-as-N |F|={neg['D_roots_as_neumann_max_|F|']:.2f}, "
        f"N-as-D resid={neg['N_roots_as_dirichlet_max_residual']:.1f} | "
        f"{controls['status']}"
    )
    footer1 = fig.text(0.5, 0.075, footer_l1, ha="center", va="bottom",
                       color=FG, fontsize=6.8, alpha=0.0)
    footer2 = fig.text(0.5, 0.022, footer_l2, ha="center", va="bottom",
                       color=FG, fontsize=6.8, alpha=0.0)

    n_intro, frames_per_pole, n_outro = 6, 6, 8
    n_frames = n_intro + len(pole_artists) * frames_per_pole + n_outro

    def update(frame: int):
        changed = []
        k = frame - n_intro
        for idx, (line, txt) in enumerate(pole_artists):
            local = k - idx * frames_per_pole
            alpha = 0.0 if local < 0 else (1.0 if local >= frames_per_pole - 1
                                           else (local + 1) / frames_per_pole)
            if line.get_alpha() != alpha:
                line.set_alpha(alpha)
                txt.set_alpha(alpha)
                changed.extend([line, txt])
        if frame >= n_frames - n_outro and footer1.get_alpha() != 1.0:
            footer1.set_alpha(1.0)
            footer2.set_alpha(1.0)
            changed.extend([footer1, footer2])
        return changed

    # blit=False: figure-level text artists (footers) only render reliably
    # in saved GIFs when each frame is drawn in full.
    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    path = _save(anim, "b6_step0_square_well_pole_reveal.gif")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Demo 2 — bridge contract chain
# ---------------------------------------------------------------------------

STAGES = [
    ("forward",
     "TransportOperator\nbound solve\n(check_forward)",
     "SSZ: null geodesics,\neikonal phase\n(own code)"),
    ("convergence",
     "resolution ladder\n(check_convergence)",
     "4N→8N rel shift\n+ overlap\ncriterion"),
    ("direction /\nparity",
     "PropagatedState\n±1 tagged\n(check_parity)",
     "Sagnac: v to -v swap;\nSSZ: NO sign-swap\n(declared rule)"),
    ("observable",
     "ObservableSet\ndirectional/scalar",
     "direction swap\nflips vectors,\npreserves scalars"),
    ("inverse",
     "check_inverse_\nround_trip",
     "Sagnac: recovers v\nfrom (t₊, t₋);\nSSZ: deferred"),
]


def demo_contract_chain() -> Path:
    cert = _load(ART / "BRIDGE_CERTIFICATE_V1.json")
    verdict = cert["verdict"]
    backends = ", ".join(sorted(cert["backends_executed"]))
    n_backends = len(cert["backends_executed"])

    fig, ax = plt.subplots(figsize=(9.8, 4.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    fig.suptitle(
        f"SSZ-Transport-Bridge — five structural contracts, one vocabulary "
        f"({n_backends} backends: {backends})",
        color=FG, fontsize=11,
    )

    boxes, captions = [], []
    for i, (name, core, note) in enumerate(STAGES):
        x0 = 0.35 + i * 1.93
        patch = FancyBboxPatch(
            (x0, 6.0), 1.62, 2.0,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.4, edgecolor=ACCENT, facecolor=ACCENT, alpha=0.0,
        )
        ax.add_patch(patch)
        label = ax.text(x0 + 0.81, 7.0, f"{name}\n\n{core}", ha="center", va="center",
                        fontsize=7.2, color=FG, alpha=0.0)
        cap = ax.text(x0 + 0.81, 5.35, note, ha="center", va="top",
                      fontsize=6.2, color=FG, alpha=0.0)
        boxes.append(patch)
        captions.append((label, cap))
        if i < len(STAGES) - 1:
            ax.annotate("", xy=(x0 + 1.98, 7.0), xytext=(x0 + 1.66, 7.0),
                        arrowprops=dict(arrowstyle="-|>", color=FG, lw=1.2))

    footer = fig.text(
        0.5, 0.12,
        f"certificate: {verdict} · no_shared_formulas="
        f"{str(cert['no_shared_formulas']).lower()} · no_skips="
        f"{str(cert['no_skips']).lower()}",
        ha="center", color=ACCENT, fontsize=8.5, alpha=0.0,
    )
    note = fig.text(
        0.5, 0.045,
        "shared TYPES and CHECK STRUCTURE only — each backend owns its physics",
        ha="center", color=FG, fontsize=7, alpha=0.0,
    )

    n_intro, frames_per_stage, n_outro = 4, 7, 8
    n_frames = n_intro + len(STAGES) * frames_per_stage + n_outro

    def update(frame: int):
        changed = []
        k = frame - n_intro
        for idx, (patch, (label, cap)) in enumerate(zip(boxes, captions, strict=True)):
            local = k - idx * frames_per_stage
            if local < 0:
                target = 0.0
            elif local >= frames_per_stage - 1:
                target = 0.28
            else:
                target = 0.28 * (local + 1) / frames_per_stage
            if patch.get_alpha() != target:
                patch.set_alpha(target)
                label.set_alpha(1.0 if target > 0.0 else 0.0)
                cap.set_alpha(1.0 if target >= 0.28 else 0.0)
                changed.extend([patch, label, cap])
        if frame >= n_frames - n_outro:
            if footer.get_alpha() != 1.0:
                footer.set_alpha(1.0)
                note.set_alpha(1.0)
                changed.extend([footer, note])
        return changed

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    path = _save(anim, "bridge_contract_chain.gif")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Demo 3 — B6 step 1: ch2 pilot pole
# ---------------------------------------------------------------------------

def demo_step1_pilot_pole() -> Path:
    pilot = _load(ART / "SSZ_B6_STEP1_CH2_PILOT_POLE_V1.json")
    status = pilot["status"]
    qnm = pilot["pole_qnm_sign"]
    scan = pilot["pole_scan_sign"]
    gates = pilot["certification_gates"]

    fig, (ax, ax_g) = plt.subplots(
        1, 2, figsize=(9.8, 4.6), gridspec_kw={"width_ratios": [1.0, 1.35]},
    )
    fig.suptitle(
        f"B6 Step 1 — ch2 pilot pole  ·  w = {qnm['re']:.5f} {qnm['im']:+.5f}i"
        f"  ·  {status}",
        color=FG, fontsize=10.5,
    )

    # left panel: complex plane near the pole
    pad = 0.03
    ax.axhline(0.0, color=GRID, lw=0.8)
    ax.plot([scan["re"]], [scan["im"]], "s", ms=9, mfc="#d29922",
            mec="none", label="pole_scan_sign (scan convention)")
    qnm_dot, = ax.plot([qnm["re"]], [qnm["im"]], "o", ms=10, mfc=ACCENT,
                       mec="none", alpha=0.0, label="pole_qnm_sign (QNM convention)")
    ax.annotate("", xy=(qnm["re"], qnm["im"]), xytext=(scan["re"], scan["im"]),
                arrowprops=dict(arrowstyle="->", color=FG, lw=1.0, alpha=0.6))
    ax.text(scan["re"] + pad, scan["im"] + pad, "scan", color="#d29922", fontsize=7)
    ax.text(qnm["re"] + pad, qnm["im"] - 2.5 * pad, "QNM sign", color=ACCENT,
            fontsize=7)
    ax.set_xlabel("Re w")
    ax.set_ylabel("Im w")
    ax.set_xlim(scan["re"] - 0.12, scan["re"] + 0.16)
    ax.set_ylim(qnm["im"] - 0.07, scan["im"] + 0.07)
    ax.grid(color=GRID, lw=0.4, alpha=0.5)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.2)

    # right panel: certification gates revealed line by line
    ax_g.axis("off")
    ax_g.set_title("certification gates", fontsize=9)
    gate_items = list(gates.items())
    n_gates = len(gate_items)
    y0, dy = 0.92, 0.92 / n_gates
    gate_texts = []
    for i, (gname, gval) in enumerate(gate_items):
        sval = str(gval)
        if len(sval) > 64:
            sval = sval[:61] + "..."
        txt = ax_g.text(
            0.02, y0 - i * dy, f"[{gname}]  {sval}",
            ha="left", va="top", color=FG, fontsize=6.4, alpha=0.0,
            family="monospace", transform=ax_g.transAxes,
        )
        gate_texts.append(txt)

    n_intro, n_mirror, frames_per_gate, n_outro = 6, 10, 4, 8
    n_frames = n_intro + n_mirror + len(gate_texts) * frames_per_gate + n_outro

    def update(frame: int):
        changed = []
        if frame >= n_intro and qnm_dot.get_alpha() != 1.0:
            ramp = min(1.0, (frame - n_intro + 1) / n_mirror)
            qnm_dot.set_alpha(ramp)
            changed.append(qnm_dot)
        k = frame - n_intro - n_mirror
        for idx, txt in enumerate(gate_texts):
            local = k - idx * frames_per_gate
            if local >= 0 and txt.get_alpha() != 1.0:
                txt.set_alpha(1.0)
                changed.append(txt)
        return changed

    anim = FuncAnimation(fig, update, frames=n_frames, blit=False)
    path = _save(anim, "b6_step1_ch2_pilot_pole.gif")
    plt.close(fig)
    return path


def main() -> int:
    _dark_theme()
    for fn in (
        demo_step0_pole_reveal,
        demo_contract_chain,
        demo_step1_pilot_pole,
    ):
        path = fn()
        print(f"wrote {path.relative_to(BRIDGE)} ({path.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
