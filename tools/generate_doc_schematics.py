from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle, Wedge


ROOT = Path(__file__).resolve().parents[1]
PHYS = ROOT / "doc" / "assets" / "figures" / "physics"
ALG = ROOT / "doc" / "assets" / "figures" / "algorithm"
DECOR = ROOT / "doc" / "assets" / "figures" / "decorative"

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
    }
)

COL = {
    "ink": "#23313d",
    "muted": "#60717f",
    "panel": "#f6f8fb",
    "blue": "#4c78a8",
    "cyan": "#72b7b2",
    "green": "#59a14f",
    "orange": "#f28e2b",
    "red": "#e15759",
    "purple": "#7b6ea8",
    "gold": "#edc948",
    "gray": "#a7b0b7",
}


def _fig(title: str, subtitle: str = "", size: tuple[float, float] = (7.4, 4.2)):
    fig, ax = plt.subplots(figsize=size)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(0.02, 0.965, title, ha="left", va="top", fontsize=10, fontweight="bold", color=COL["ink"])
    if subtitle:
        ax.text(0.02, 0.91, subtitle, ha="left", va="top", fontsize=7, color=COL["muted"])
    return fig, ax


def _bg(ax, hue: str = "blue"):
    base = np.linspace(0, 1, 320)
    grad = np.outer(np.ones_like(base), base)
    palettes = {
        "blue": ((0.93, 0.97, 1.00), (0.99, 0.99, 0.96)),
        "purple": ((0.96, 0.94, 1.00), (0.98, 0.99, 0.95)),
        "orange": ((1.00, 0.96, 0.90), (0.94, 0.98, 1.00)),
        "green": ((0.92, 0.98, 0.94), (0.99, 0.97, 0.92)),
    }
    c0, c1 = palettes[hue]
    img = np.zeros((320, 320, 3))
    for k in range(3):
        img[:, :, k] = c0[k] * (1 - grad) + c1[k] * grad
    ax.imshow(img, extent=(0, 1, 0, 1), origin="lower", zorder=-10, aspect="auto")


def _box(ax, x, y, w, h, text, fc="panel", ec="ink", color=None, fs=7, alpha=0.96):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.022",
        fc=COL.get(fc, fc),
        ec=COL.get(ec, ec),
        lw=0.9,
        alpha=alpha,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color or COL["ink"], linespacing=1.15)
    return patch


def _arrow(ax, a, b, color="ink", rad=0.0, lw=1.2, ms=10, alpha=1.0):
    arr = FancyArrowPatch(
        a,
        b,
        arrowstyle="-|>",
        mutation_scale=ms,
        lw=lw,
        color=COL.get(color, color),
        connectionstyle=f"arc3,rad={rad}",
        alpha=alpha,
    )
    ax.add_patch(arr)
    return arr


def _label(ax, x, y, text, fs=7, color="ink", ha="center", va="center", weight=None):
    ax.text(x, y, text, fontsize=fs, color=COL.get(color, color), ha=ha, va=va, fontweight=weight)


def _particles(ax, center=(0.5, 0.5), radius=0.35, n=80, color="blue", seed=0):
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2 * np.pi, n)
    rr = radius * np.sqrt(rng.uniform(0, 1, n))
    x = center[0] + rr * np.cos(theta)
    y = center[1] + rr * np.sin(theta)
    ax.scatter(x, y, s=rng.uniform(4, 14, n), c=COL[color], alpha=0.12, lw=0)


def _stars(ax, n=180, seed=0):
    rng = np.random.default_rng(seed)
    ax.scatter(
        rng.uniform(0, 1, n),
        rng.uniform(0, 1, n),
        s=rng.uniform(2, 14, n),
        c="white",
        alpha=rng.uniform(0.08, 0.35, n),
        lw=0,
        zorder=-2,
    )


def _save(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf"):
        target = stem.with_suffix(f".{ext}")
        if target.exists():
            target.unlink()
        fig.savefig(target, bbox_inches="tight", facecolor="white")
    png = stem.with_suffix(".png")
    if png.exists():
        png.unlink()
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    tiff = stem.with_suffix(".tiff")
    if tiff.exists():
        tiff.unlink()
    fig.savefig(
        tiff,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


def physical_chain():
    fig, ax = _fig("ASGARD physical chain", "Local shock physics is radius-ordered; observer quantities are projected last.")
    _bg(ax, "blue")
    _particles(ax, (0.52, 0.50), 0.44, 120, "blue", 1)
    xs = [0.05, 0.23, 0.41, 0.59, 0.77]
    texts = ["relativistic\nejecta", "FS / RS\nshocks", "e$^-$ / p\ntransport", "photons +\nsecondaries", "EATS\nprojection"]
    cols = ["purple", "orange", "green", "blue", "cyan"]
    for x, t, c in zip(xs, texts, cols):
        _box(ax, x, 0.49, 0.13, 0.15, t, fc=c, ec=c, color="white")
    for a, b in zip(xs[:-1], xs[1:]):
        _arrow(ax, (a + 0.13, 0.565), (b, 0.565), "ink")
    _label(ax, 0.51, 0.29, r"$R_i\rightarrow\{\Gamma_i,N_e,N_p,n_\gamma\}\rightarrow F_\nu(t_{\rm obs})$", 9)
    _save(fig, PHYS / "physical_chain")


def decorative_physics_header():
    fig, ax = _fig(
        "Afterglow physical geometry",
        "Central engine, jet, reverse shock and forward shock all propagate left-to-right.",
        size=(9.0, 3.6),
    )
    # Dark editorial background.
    x = np.linspace(0, 1, 420)
    y = np.linspace(0, 1, 240)
    xx, yy = np.meshgrid(x, y)
    img = np.zeros((240, 420, 3))
    img[..., 0] = 0.025 + 0.07 * xx + 0.05 * np.exp(-((xx - 0.75) ** 2 + (yy - 0.50) ** 2) / 0.05)
    img[..., 1] = 0.050 + 0.10 * xx + 0.10 * np.exp(-((xx - 0.72) ** 2 + (yy - 0.50) ** 2) / 0.06)
    img[..., 2] = 0.100 + 0.18 * xx + 0.18 * np.exp(-((xx - 0.70) ** 2 + (yy - 0.50) ** 2) / 0.08)
    ax.imshow(img, extent=(0, 1, 0, 1), origin="lower", aspect="auto", zorder=-10)
    _stars(ax, 160, 42)

    # Engine and right-going jet cone.
    engine = (0.12, 0.50)
    ax.add_patch(Circle(engine, 0.045, fc="#fff2c6", ec="#fff9e5", lw=1.0, alpha=0.95, zorder=4))
    ax.add_patch(Circle(engine, 0.085, fc="#edc948", ec="none", alpha=0.10, zorder=2))
    upper = np.array([[0.15, 0.52], [0.73, 0.68], [0.82, 0.58], [0.20, 0.53]])
    lower = np.array([[0.15, 0.48], [0.73, 0.32], [0.82, 0.42], [0.20, 0.47]])
    ax.add_patch(Polygon(upper, closed=True, fc="#72b7b2", ec="none", alpha=0.18, zorder=1))
    ax.add_patch(Polygon(lower, closed=True, fc="#72b7b2", ec="none", alpha=0.18, zorder=1))
    ax.plot([0.13, 0.86], [0.50, 0.50], color="#b8fff4", lw=1.0, alpha=0.55, zorder=2)

    # Reverse shock shell behind the leading shock, both convex and moving right.
    ax.add_patch(Wedge((0.70, 0.50), 0.31, -38, 38, width=0.040, fc="#f28e2b", ec="#ffd19b", lw=1.0, alpha=0.78, zorder=3))
    ax.add_patch(Wedge((0.80, 0.50), 0.36, -34, 34, width=0.050, fc="#72d6ff", ec="#d9f6ff", lw=1.2, alpha=0.78, zorder=4))
    ax.add_patch(Wedge((0.84, 0.50), 0.43, -30, 30, width=0.020, fc="#dff8ff", ec="none", alpha=0.25, zorder=2))

    # Particle/photon trails strictly left-to-right.
    rng = np.random.default_rng(7)
    for _ in range(42):
        y0 = 0.50 + rng.normal(0, 0.10)
        x0 = rng.uniform(0.19, 0.70)
        x1 = x0 + rng.uniform(0.08, 0.18)
        y1 = y0 + rng.normal(0, 0.015)
        ax.plot([x0, x1], [y0, y1], color="#e7f7ff", lw=rng.uniform(0.3, 0.8), alpha=rng.uniform(0.12, 0.45), zorder=2)
    for x0 in [0.26, 0.44, 0.62]:
        _arrow(ax, (x0, 0.20), (x0 + 0.15, 0.20), color="#e7f7ff", ms=9, alpha=0.75)
    _label(ax, 0.12, 0.22, "engine", 7, "#f7f1d2")
    _label(ax, 0.50, 0.22, "jet flow", 7, "#d8f7ff")
    _label(ax, 0.72, 0.78, "reverse shock", 7, "#ffd19b")
    _label(ax, 0.84, 0.82, "forward shock", 7, "#d9f6ff")
    _save(fig, DECOR / "physics_afterglow_header")


def decorative_algorithm_header():
    fig, ax = _fig(
        "Algorithmic data flow",
        "Concrete state objects, transport grids and observer projection are separated left-to-right.",
        size=(9.0, 3.6),
    )
    _bg(ax, "blue")
    # Configuration cards.
    for i, (x, y, c) in enumerate([(0.06, 0.66, "blue"), (0.08, 0.55, "green"), (0.10, 0.44, "purple")]):
        ax.add_patch(Rectangle((x, y), 0.12, 0.075, fc=COL[c], ec="white", lw=0.8, alpha=0.86, zorder=2))
    # Radius shells.
    for r in [0.08, 0.12, 0.16]:
        ax.add_patch(Wedge((0.30, 0.53), r, -55, 55, width=0.012, fc="none", ec=COL["cyan"], lw=1.4, alpha=0.85))
    # Electron finite volume mesh.
    for i in range(6):
        for j in range(4):
            color = "#e5f1ec" if (i + j) % 2 == 0 else "#d2e8df"
            ax.add_patch(Rectangle((0.43 + i * 0.035, 0.39 + j * 0.045), 0.032, 0.040, fc=color, ec="white", lw=0.6))
    for offset in [0.00, 0.05, 0.10]:
        _arrow(ax, (0.45 + offset, 0.60), (0.55 + offset, 0.43), color="green", ms=7, alpha=0.85)
    # Solver and photon/hadronic panels.
    _box(ax, 0.68, 0.55, 0.10, 0.10, "", "orange", "orange")
    for k in range(3):
        ax.plot([0.695, 0.765], [0.575 + 0.018 * k, 0.575 + 0.018 * k], color="white", lw=1.0)
    spectrum_x = np.linspace(0.67, 0.80, 80)
    ax.plot(spectrum_x, 0.38 + 0.10 * np.exp(-((spectrum_x - 0.73) / 0.025) ** 2), color=COL["blue"], lw=1.6)
    ax.add_patch(Rectangle((0.65, 0.32), 0.18, 0.18, fc="none", ec=COL["blue"], lw=1.0, alpha=0.65))
    # Observer grid.
    for i in range(5):
        ax.plot([0.88, 0.96], [0.34 + i * 0.045, 0.34 + i * 0.045], color=COL["purple"], lw=0.8)
        ax.plot([0.88 + i * 0.02, 0.88 + i * 0.02], [0.34, 0.52], color=COL["purple"], lw=0.8)
    for a, b in [((0.20, 0.58), (0.27, 0.55)), ((0.36, 0.53), (0.43, 0.50)), ((0.64, 0.51), (0.68, 0.60)), ((0.78, 0.43), (0.88, 0.43))]:
        _arrow(ax, a, b, color="ink", ms=9)
    _label(ax, 0.12, 0.31, "configuration", 7)
    _label(ax, 0.30, 0.31, "radius shells", 7)
    _label(ax, 0.52, 0.31, "transport grid", 7)
    _label(ax, 0.73, 0.27, "local kernels", 7)
    _label(ax, 0.92, 0.27, "observer grid", 7)
    _save(fig, DECOR / "algorithm_flow_header")


def spacetime():
    fig, ax = _fig("Coordinates and observer mapping", "Transport uses R; Doppler and equal-arrival-time geometry set the observed light curve.")
    _bg(ax, "purple")
    ax.plot([0.10, 0.88], [0.18, 0.82], lw=2.0, color=COL["blue"])
    for i, r in enumerate(np.linspace(0.2, 0.78, 5)):
        y = 0.10 + 0.82 * r
        ax.add_patch(Circle((r, y), 0.022, fc=COL["orange"], ec="white", lw=0.6))
        _arrow(ax, (r, y), (r + 0.08, y + 0.04), "orange", ms=8)
    _label(ax, 0.29, 0.78, r"$dt^\prime/dR=(\beta\Gamma c)^{-1}$", 9)
    _label(ax, 0.67, 0.31, r"$t_{\rm obs}=(1+z)(t_{\rm lab}-R\mu/c)$", 9)
    _label(ax, 0.64, 0.70, r"$\delta=[\Gamma(1-\beta\mu)]^{-1}$", 9)
    _save(fig, PHYS / "spacetime_observer_mapping")


def medium_jet():
    fig, ax = _fig("Medium and jet structure", "Density profile and angular energy profile define each blast-wave patch.")
    _bg(ax, "green")
    rr = np.linspace(0, 1, 240)
    ax.plot(0.06 + 0.38 * rr, 0.74 + 0 * rr, color=COL["blue"], lw=2)
    ax.plot(0.06 + 0.38 * rr, 0.30 + 0.42 / (1 + 8 * rr), color=COL["orange"], lw=2)
    ax.plot(0.06 + 0.38 * rr, 0.28 + 0.22 * np.exp(-((rr - 0.58) / 0.08) ** 2), color=COL["green"], lw=2)
    _label(ax, 0.24, 0.82, r"$n_{\rm ISM}=n_0$", color="blue")
    _label(ax, 0.24, 0.53, r"$n_w=3\times10^{35}A_*R^{-2}$", color="orange")
    theta = np.linspace(-0.55, 0.55, 180)
    e = np.exp(-0.5 * (theta / 0.18) ** 2)
    x = 0.62 + 0.29 * theta / 0.55
    y = 0.23 + 0.52 * e
    ax.fill_between(x, 0.23, y, color=COL["purple"], alpha=0.16)
    ax.plot(x, y, color=COL["purple"], lw=2.2)
    _label(ax, 0.76, 0.82, r"$E_{\rm iso}(\theta),\,\Gamma_0(\theta)$")
    _save(fig, PHYS / "medium_jet_structure")


def forward_dynamics():
    fig, ax = _fig("Forward-shock dynamics", "Swept-up mass converts initial ejecta energy into decelerating blast-wave motion.")
    _bg(ax, "orange")
    _box(ax, 0.07, 0.59, 0.18, 0.13, r"$n(R)$", "blue", "blue", "white")
    _box(ax, 0.07, 0.34, 0.18, 0.13, r"$E_{\rm iso},\,\Gamma_0$", "purple", "purple", "white")
    _box(ax, 0.36, 0.45, 0.27, 0.16, r"$dM/dR=4\pi R^2nm_p$" + "\n" + r"$E\simeq C_k\Gamma^2Mc^2$")
    _box(ax, 0.75, 0.47, 0.16, 0.13, r"$\Gamma(R)$" + "\n" + r"$t_{\rm obs}(R)$", "orange", "orange", "white")
    _arrow(ax, (0.25, 0.65), (0.36, 0.55))
    _arrow(ax, (0.25, 0.40), (0.36, 0.51))
    _arrow(ax, (0.63, 0.53), (0.75, 0.53))
    _label(ax, 0.50, 0.24, r"$\Gamma\propto R^{-(3-k)/2}$,  $\Gamma\propto t_{\rm obs}^{-(3-k)/2(4-k)}$", 8)
    _save(fig, PHYS / "forward_shock_dynamics")


def electron_transport():
    fig, ax = _fig("Electron injection and transport", "Injection, cooling and conservative advection determine the electron spectrum.")
    _bg(ax, "green")
    _box(ax, 0.06, 0.56, 0.21, 0.15, r"$Q_e\propto\gamma^{-p}e^{-\gamma/\gamma_{\max}}$", "green", "green", "white")
    _box(ax, 0.39, 0.56, 0.20, 0.15, r"$\dot\gamma_R$" + "\n" + "ad + syn + IC", "orange", "orange", "white")
    _box(ax, 0.72, 0.56, 0.20, 0.15, r"$N_e(\gamma,R)$", "blue", "blue", "white")
    _arrow(ax, (0.27, 0.635), (0.39, 0.635))
    _arrow(ax, (0.59, 0.635), (0.72, 0.635))
    _label(ax, 0.50, 0.34, r"$\partial_RN_e+\partial_\gamma(\dot\gamma_RN_e)=Q_{e,R}$", 10)
    _label(ax, 0.50, 0.22, r"$\dot\gamma_{\rm syn}^\prime=-\sigma_TB^{\prime2}\gamma^2/(6\pi m_ec)$", 8)
    _save(fig, PHYS / "electron_transport")


def synch_ssc():
    fig, ax = _fig("Synchrotron, SSA and SSC", "One electron distribution produces local photons, absorption and inverse-Compton emission.")
    _bg(ax, "blue")
    _box(ax, 0.06, 0.50, 0.16, 0.15, r"$N_e(\gamma)$", "green", "green", "white")
    _box(ax, 0.31, 0.64, 0.20, 0.12, "synchrotron\n" + r"$L_\nu^\prime$", "blue", "blue", "white")
    _box(ax, 0.31, 0.38, 0.20, 0.12, "SSA\n" + r"$\tau_\nu=\alpha_\nu\ell$", "orange", "orange", "white")
    _box(ax, 0.64, 0.51, 0.21, 0.14, "SSC / IC\n" + r"$e^-+\gamma$", "purple", "purple", "white")
    _arrow(ax, (0.22, 0.58), (0.31, 0.70))
    _arrow(ax, (0.22, 0.55), (0.31, 0.44))
    _arrow(ax, (0.51, 0.70), (0.64, 0.59))
    _arrow(ax, (0.51, 0.44), (0.64, 0.56))
    _label(ax, 0.50, 0.23, r"$L_\nu^\prime=\int N_eP_\nu^\prime d\gamma$,  $S_\nu=(1-e^{-\tau_\nu})/\tau_\nu$", 8)
    _save(fig, PHYS / "synch_ssa_ssc")


def hadronic_transport():
    fig, ax = _fig("Hadronic transport", "Proton losses and secondary production share the same local photon and baryon targets.")
    _bg(ax, "purple")
    _box(ax, 0.06, 0.48, 0.16, 0.16, r"$Q_p\rightarrow N_p$", "purple", "purple", "white")
    for text, color, xy in [
        ("p syn", "blue", (0.36, 0.74)),
        (r"p$\gamma$", "orange", (0.36, 0.55)),
        ("BH pairs", "green", (0.36, 0.36)),
        ("pp", "red", (0.66, 0.63)),
        ("neutrino", "cyan", (0.66, 0.43)),
    ]:
        _box(ax, xy[0], xy[1], 0.18, 0.10, text, color, color, "white")
        _arrow(ax, (0.22, 0.56), (xy[0], xy[1] + 0.05), rad=0.04)
    _label(ax, 0.50, 0.18, r"$t_{p\gamma}^{-1}\sim c\int\sigma_{p\gamma}\kappa\,n_\gamma\,d\epsilon$", 9)
    _save(fig, PHYS / "hadronic_transport")


def joint_feedback():
    fig, ax = _fig("Joint secondary feedback", "Electrons, protons and photons close on the same R grid.")
    _bg(ax, "orange")
    _box(ax, 0.13, 0.57, 0.18, 0.12, r"$N_e$", "green", "green", "white")
    _box(ax, 0.41, 0.68, 0.18, 0.12, r"$n_\gamma$", "blue", "blue", "white")
    _box(ax, 0.70, 0.57, 0.18, 0.12, r"$N_p$", "purple", "purple", "white")
    _box(ax, 0.38, 0.31, 0.24, 0.12, r"$Q_{e,\rm sec}$ + photon sink", "orange", "orange", "white")
    _arrow(ax, (0.31, 0.63), (0.41, 0.73))
    _arrow(ax, (0.59, 0.74), (0.70, 0.64))
    _arrow(ax, (0.79, 0.57), (0.62, 0.39))
    _arrow(ax, (0.38, 0.38), (0.22, 0.57))
    _arrow(ax, (0.50, 0.43), (0.50, 0.68), rad=-0.20)
    _label(ax, 0.50, 0.18, r"$N_e=E[Q_{\rm prim}+Q_{\rm sec},n_\gamma],\quad N_p=H[Q_p,n_\gamma]$", 8)
    _save(fig, PHYS / "joint_feedback")


def pair_cascade():
    fig, ax = _fig("Gamma-gamma pair cascade", "Absorbed gamma rays create pairs whose synchrotron seed is carried shell by shell.")
    _bg(ax, "purple")
    _box(ax, 0.06, 0.61, 0.18, 0.11, r"$\gamma_{\rm HE}$", "purple", "purple", "white")
    _box(ax, 0.06, 0.37, 0.18, 0.11, r"$\gamma_{\rm target}$", "blue", "blue", "white")
    _box(ax, 0.38, 0.49, 0.19, 0.13, r"$\gamma\gamma\rightarrow e^+e^-$", "orange", "orange", "white")
    _box(ax, 0.69, 0.49, 0.20, 0.13, "pair synch\nnext shell", "green", "green", "white")
    _arrow(ax, (0.24, 0.665), (0.38, 0.56))
    _arrow(ax, (0.24, 0.425), (0.38, 0.54))
    _arrow(ax, (0.57, 0.555), (0.69, 0.555))
    _label(ax, 0.49, 0.25, r"$E_1E_2(1-\cos\psi)\geq2(m_ec^2)^2$", 9)
    _save(fig, PHYS / "pair_cascade")


def reverse_shock():
    fig, ax = _fig("Reverse shock and density-bump secondary RS", "Region 3 emission uses newly dissipated energy rather than recycled FS heat.")
    _bg(ax, "orange")
    for x, text, color in [
        (0.08, "1\ncold\nmedium", "gray"),
        (0.30, "2\ntransmitted\nFS", "blue"),
        (0.52, "3\nsecondary\nRS", "orange"),
        (0.74, "4\nhot old\nshell", "purple"),
    ]:
        _box(ax, x, 0.45, 0.15, 0.22, text, color, color, "white")
    for x in (0.25, 0.49, 0.71):
        _arrow(ax, (x, 0.56), (x + 0.04, 0.56), ms=8)
    _label(ax, 0.50, 0.32, r"$\Gamma_2=\Gamma_3=\Gamma_c,\quad p_2(\Gamma_c)=p_3(\Gamma_c)$", 9)
    _label(ax, 0.50, 0.20, r"$\gamma_{m,3}=1+\frac{\epsilon_{e,3}}{\xi_{N,3}}\frac{p_3-2}{p_3-1}\frac{u_{\rm diss,3}}{n_3m_ec^2}$", 8)
    _save(fig, PHYS / "reverse_shock_secondary")


def chi_projection():
    fig, ax = _fig(r"$\chi$-resolved finite-thickness projection", "Only FS synchrotron + SSA uses this finite-thickness EATS contract.")
    _bg(ax, "blue")
    center = (0.35, 0.44)
    for rad, alpha in [(0.33, 0.13), (0.25, 0.20), (0.17, 0.27)]:
        ax.add_patch(Wedge(center, rad, 18, 80, width=0.055, fc=COL["blue"], alpha=alpha, ec=COL["blue"], lw=0.8))
    _arrow(ax, (0.55, 0.58), (0.83, 0.72), "orange")
    _label(ax, 0.30, 0.78, r"$\chi=1+2(4-k)\Gamma_{\rm sh}^2(1-r/R)$", 8)
    _label(ax, 0.76, 0.77, "observer")
    _label(ax, 0.68, 0.35, r"$F_\nu=\sum W_R W_\Omega W_\chi\delta^3L_{\nu'}^\prime S_\nu$", 8)
    _save(fig, PHYS / "chi_eats_projection")


def alg_pipeline():
    fig, ax = _fig("Public API to solve state", "User semantics are normalized before numerical kernels run.")
    _bg(ax, "blue")
    xs = [0.04, 0.20, 0.38, 0.56, 0.74]
    texts = ["Model\nAPI", "Runtime\nConfig", "Simulation\nSetup", "SolveState", "Flux /\ndetails"]
    for x, t, c in zip(xs, texts, ["blue", "cyan", "green", "orange", "purple"]):
        _box(ax, x, 0.50, 0.13, 0.15, t, c, c, "white")
    for a, b in zip(xs[:-1], xs[1:]):
        _arrow(ax, (a + 0.13, 0.575), (b, 0.575))
    _label(ax, 0.50, 0.30, "projection reads transport state; it never rewrites particle spectra", 8)
    _save(fig, ALG / "api_to_solve_state")


def alg_grid():
    fig, ax = _fig("Log-grid conservative variables", "Energy and frequency coordinates carry explicit Jacobians.")
    _bg(ax, "green")
    xs = np.linspace(0.12, 0.88, 9)
    for x in xs:
        ax.plot([x, x], [0.45, 0.62], color=COL["gray"], lw=1)
    ax.plot([xs[0], xs[-1]], [0.535, 0.535], color=COL["ink"], lw=1.2)
    _label(ax, 0.50, 0.71, r"$x=\log_{10}\gamma,\quad N_x=dN/dx=\gamma\ln10\,dN/d\gamma$", 9)
    _label(ax, 0.50, 0.30, r"$L_y=dL/dy=\nu\ln10\,L_\nu,\quad y=\log_{10}\nu$", 9)
    _save(fig, ALG / "log_grid_jacobian")


def alg_rk():
    fig, ax = _fig("Dynamics RK update with event splitting", "Shock crossings or density events split the step before integration.")
    _bg(ax, "orange")
    ax.plot([0.08, 0.90], [0.55, 0.55], color=COL["ink"], lw=1.5)
    for x, t, c in [(0.12, r"$q_n$", "ink"), (0.48, r"$q_*$", "red"), (0.86, r"$q_{n+1}$", "ink")]:
        ax.plot([x, x], [0.50, 0.62], color=COL[c], lw=1.5)
        _label(ax, x, 0.45, t)
    _arrow(ax, (0.14, 0.68), (0.46, 0.68), "blue")
    _arrow(ax, (0.50, 0.68), (0.84, 0.68), "orange")
    _label(ax, 0.50, 0.27, r"$y_{n+1}=y_n+\frac{h}{6}(k_1+2k_2+2k_3+k_4)$, on each physical branch", 8)
    _save(fig, ALG / "rk_event_splitting")


def alg_fullhide():
    fig, ax = _fig("fullhide finite-volume electron update", "Implicit upwind transport solves a linear system, not a smoothed output.")
    _bg(ax, "purple")
    for i in range(6):
        for j in range(4):
            ax.add_patch(Rectangle((0.12 + i * 0.09, 0.30 + j * 0.09), 0.083, 0.083, fc="#eef3f7", ec="white", lw=0.8))
    _arrow(ax, (0.62, 0.58), (0.80, 0.58), "orange")
    _box(ax, 0.76, 0.48, 0.17, 0.16, r"$A N^{n+1}=b$" + "\ntridiagonal", "orange", "orange", "white")
    _label(ax, 0.33, 0.75, r"$N_i^{n+1}=N_i^n-\Delta R(F_{i+1/2}-F_{i-1/2})/\Delta x+\Delta RQ_i$", 8)
    _label(ax, 0.33, 0.18, r"cooling flow in log-$\gamma$ cells", 8)
    _save(fig, ALG / "fullhide_finite_volume")


def alg_substep():
    fig, ax = _fig("Substeps and adaptive transport error", "Accuracy is judged on transport state, not post-processed flux smoothness.")
    _bg(ax, "blue")
    ax.plot([0.12, 0.88], [0.62, 0.62], color=COL["ink"], lw=1.2)
    for x in np.linspace(0.12, 0.88, 7):
        ax.plot([x, x], [0.58, 0.66], color=COL["ink"], lw=1)
    _arrow(ax, (0.12, 0.42), (0.88, 0.42), "blue")
    _arrow(ax, (0.12, 0.30), (0.50, 0.30), "orange")
    _arrow(ax, (0.50, 0.30), (0.88, 0.30), "orange")
    _label(ax, 0.50, 0.73, r"$\delta R=(R_{n+1}-R_n)/N_{\rm sub}$", 9)
    _label(ax, 0.50, 0.18, r"$\epsilon=\|N_{1/2,1/2}-N_1\|_w/\|N_{1/2,1/2}\|_w$", 9)
    _save(fig, ALG / "substep_error_control")


def alg_photon():
    fig, ax = _fig("PhotonFieldState construction", "Local target photons are not observer luminosities.")
    _bg(ax, "green")
    _box(ax, 0.06, 0.52, 0.18, 0.14, r"$L_\nu^\prime$" + "\nsynch seed", "blue", "blue", "white")
    for y, text, c in [(0.70, "target\n" + r"$n_\nu^\prime$", "green"), (0.50, "absorption\n" + r"$\tau_\nu,S_\nu$", "orange"), (0.30, "SSC seed", "purple")]:
        _box(ax, 0.38, y, 0.20, 0.12, text, c, c, "white")
        _arrow(ax, (0.24, 0.59), (0.38, y + 0.06))
    _box(ax, 0.72, 0.50, 0.20, 0.14, "observer\n" + r"$F_\nu$", "gray", "gray", "white")
    _arrow(ax, (0.58, 0.56), (0.72, 0.57))
    _label(ax, 0.50, 0.18, r"$n_\nu^\prime\simeq L_\nu^\prime t_{\rm esc}^\prime/(h\nu^\prime V^\prime)$ needs geometry", 8)
    _save(fig, ALG / "photon_field_state")


def alg_ssc():
    fig, ax = _fig("SSC spectrum and IC cooling order", "Cooling and emitted SSC photons must use the same seed field.")
    _bg(ax, "orange")
    _box(ax, 0.08, 0.58, 0.18, 0.12, r"$N_e$", "green", "green", "white")
    _box(ax, 0.38, 0.70, 0.20, 0.12, "IC cooling\n" + r"$\dot\gamma_{\rm IC}$", "orange", "orange", "white")
    _box(ax, 0.38, 0.42, 0.20, 0.12, "SSC photons\n" + r"$j_{\nu_s}$", "purple", "purple", "white")
    _box(ax, 0.70, 0.56, 0.20, 0.12, "energy budget\ncheck", "blue", "blue", "white")
    _label(ax, 0.25, 0.33, r"$n_\nu^\prime$ seed", 9, "blue")
    _arrow(ax, (0.26, 0.64), (0.38, 0.76))
    _arrow(ax, (0.26, 0.62), (0.38, 0.48))
    _arrow(ax, (0.58, 0.76), (0.70, 0.63))
    _arrow(ax, (0.58, 0.48), (0.70, 0.60))
    _save(fig, ALG / "ssc_ic_order")


def alg_hadronic():
    fig, ax = _fig("Hadronic transport algorithm", "Microphysics rates are converted to R-coordinate shell sources with Jacobians.")
    _bg(ax, "purple")
    _box(ax, 0.05, 0.56, 0.18, 0.13, "photon target\n" + r"$n_\gamma$", "blue", "blue", "white")
    _box(ax, 0.05, 0.34, 0.18, 0.13, "proton grid\n" + r"$N_p$", "purple", "purple", "white")
    _box(ax, 0.38, 0.45, 0.22, 0.16, "Fortran\nhadronic kernels", "orange", "orange", "white")
    _box(ax, 0.73, 0.45, 0.18, 0.16, "losses, pairs,\nphotons, neutrinos", "green", "green", "white")
    _arrow(ax, (0.23, 0.62), (0.38, 0.55))
    _arrow(ax, (0.23, 0.40), (0.38, 0.50))
    _arrow(ax, (0.60, 0.53), (0.73, 0.53))
    _label(ax, 0.50, 0.22, r"$Q_{x,R}=\gamma\ln10\,\dot N_\gamma^\prime\,dt^\prime/dR$", 9)
    _save(fig, ALG / "hadronic_algorithm")


def alg_joint():
    fig, ax = _fig("Joint feedback fixed-point iteration", "Each shell closes electrons, protons and photons before observer projection.")
    _bg(ax, "orange")
    centers = [(0.50, 0.77), (0.76, 0.54), (0.50, 0.31), (0.24, 0.54)]
    texts = [r"$n_\gamma^{(m)}$", r"$N_p^{(m+1)}$", r"$Q_{e,\rm sec}^{(m+1)}$", r"$N_e^{(m+1)}$"]
    cols = ["blue", "purple", "orange", "green"]
    for (x, y), t, c in zip(centers, texts, cols):
        _box(ax, x - 0.085, y - 0.045, 0.17, 0.09, t, c, c, "white")
    for a, b in zip(centers, centers[1:] + centers[:1]):
        _arrow(ax, a, b, rad=0.16, ms=8)
    _label(ax, 0.50, 0.18, "paired source and photon sink come from the same operator", 9)
    _save(fig, ALG / "joint_feedback_iteration")


def alg_pair():
    fig, ax = _fig("Pair-cascade shell sequence", "Iteration follows shell order; IC-mediated cascade is outside the current contract.")
    _bg(ax, "blue")
    for i, x in enumerate([0.10, 0.30, 0.50, 0.70]):
        _box(ax, x, 0.50, 0.14, 0.14, f"shell {i+1}\n" + r"$\gamma\gamma$ + pair syn", "blue" if i % 2 == 0 else "cyan", "blue", "white", fs=6)
        if i < 3:
            _arrow(ax, (x + 0.14, 0.57), (x + 0.20, 0.57))
    _label(ax, 0.50, 0.29, r"$n_{\gamma,i}^{new}=n_{\gamma,i}^{surv}+n_{{pair\,syn},i-1}$", 9)
    _save(fig, ALG / "pair_cascade_iteration")


def alg_eats():
    fig, ax = _fig("EATS interpolation algorithm", "Shell luminosity is Doppler-shifted and deposited on observer-time cells.")
    _bg(ax, "purple")
    for r in [0.18, 0.28, 0.38]:
        ax.add_patch(Wedge((0.30, 0.35), r, 18, 80, width=0.015, ec=COL["blue"], fc="none", lw=1.5))
    for p in [(0.38, 0.47), (0.42, 0.56), (0.45, 0.64)]:
        _arrow(ax, p, (0.80, 0.76), "orange", ms=8)
    _box(ax, 0.66, 0.20, 0.25, 0.12, "time-frequency\nobserver grid", "green", "green", "white")
    _label(ax, 0.52, 0.83, r"$t_{\rm obs}=(1+z)(t_{\rm lab}-R\mu/c)$", 9)
    _label(ax, 0.52, 0.12, r"$F_\nu(T_a)=\sum W I_a(t_{\rm obs})\delta^3L_{\nu^\prime}^\prime$", 8)
    _save(fig, ALG / "eats_interpolation")


def alg_chi():
    fig, ax = _fig("chi-resolved projection algorithm", "Transport chi cells are conservatively remapped to shell-local projection cells.")
    _bg(ax, "green")
    for i, x in enumerate(np.linspace(0.10, 0.42, 5)):
        ax.add_patch(Rectangle((x, 0.54), 0.055, 0.16, fc=COL["blue"], alpha=0.18 + 0.12 * i, ec="white"))
    for i, x in enumerate(np.linspace(0.58, 0.86, 4)):
        ax.add_patch(Rectangle((x, 0.38), 0.06, 0.24, fc=COL["orange"], alpha=0.25 + 0.10 * i, ec="white"))
    _arrow(ax, (0.45, 0.60), (0.58, 0.52))
    _label(ax, 0.26, 0.76, "transport chi grid")
    _label(ax, 0.74, 0.66, "projection chi grid")
    _label(ax, 0.50, 0.23, r"conserve $\sum P\Delta\chi$ and cell-integrated optical depth", 9)
    _save(fig, ALG / "chi_projection_algorithm")


def alg_patch_cache():
    fig, ax = _fig("Structured patches, cache and validation", "Patch aggregation is separate from cold solve versus warm query timing.")
    _bg(ax, "orange")
    for x, y in [(0.12, 0.58), (0.24, 0.68), (0.24, 0.48), (0.36, 0.58)]:
        ax.add_patch(Circle((x, y), 0.045, fc=COL["cyan"], ec="white", lw=1))
    _box(ax, 0.50, 0.52, 0.18, 0.15, r"$\sum_jF_{\nu,j}\Delta\Omega_j$", "purple", "purple", "white")
    _box(ax, 0.76, 0.62, 0.16, 0.11, "cold solve", "orange", "orange", "white")
    _box(ax, 0.76, 0.38, 0.16, 0.11, "warm query", "green", "green", "white")
    _arrow(ax, (0.40, 0.58), (0.50, 0.60))
    _arrow(ax, (0.68, 0.60), (0.76, 0.67))
    _arrow(ax, (0.68, 0.57), (0.76, 0.43))
    _label(ax, 0.50, 0.22, r"$t_{\rm cold}=t_{\rm dyn}+t_e+t_\gamma+t_p+t_{\rm proj}$;  $t_{\rm warm}\simeq t_{\rm query/proj}$", 8)
    _save(fig, ALG / "patch_cache_validation")


def main():
    funcs = [
        physical_chain,
        spacetime,
        medium_jet,
        forward_dynamics,
        electron_transport,
        synch_ssc,
        hadronic_transport,
        joint_feedback,
        pair_cascade,
        reverse_shock,
        chi_projection,
        alg_pipeline,
        alg_grid,
        alg_rk,
        alg_fullhide,
        alg_substep,
        alg_photon,
        alg_ssc,
        alg_hadronic,
        alg_joint,
        alg_pair,
        alg_eats,
        alg_chi,
        alg_patch_cache,
    ]
    for func in funcs:
        func()
    print(f"generated {len(list(PHYS.glob('*.png')))} physics PNGs and {len(list(ALG.glob('*.png')))} algorithm PNGs")


if __name__ == "__main__":
    main()
