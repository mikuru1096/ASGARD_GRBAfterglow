from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Wedge
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper"
FIG_DIR = PAPER / "figures"
DATA_DIR = PAPER / "source_data"
BENCH = DATA_DIR / "benchmarks"

BLUE = "#24588D"
TEAL = "#2E8B8B"
VIOLET = "#7750A6"
RED = "#B34A48"
GOLD = "#B8860B"
GRAY = "#4A4A4A"
LIGHT = "#E6E6E6"
PALE = "#F5F7F9"
COLORS = (BLUE, TEAL, VIOLET, RED, GOLD, GRAY)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def values(table: list[dict[str, str]], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in table], dtype=float)


def save_pub(fig: plt.Figure, stem: str) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    svg_path = FIG_DIR / f"{stem}.svg"
    fig.savefig(svg_path)
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
    fig.savefig(FIG_DIR / f"{stem}.pdf")
    fig.savefig(FIG_DIR / f"{stem}.tiff", dpi=600)
    plt.close(fig)


def panel(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.04, label, transform=ax.transAxes, fontsize=9, fontweight="bold")


def finish(ax: plt.Axes, xlabel: str = "", ylabel: str = "", title: str = "") -> None:
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(color=LIGHT, lw=0.5, which="both", zorder=0)


def flow_box(ax: plt.Axes, x: float, text: str, color: str) -> None:
    ax.add_patch(FancyBboxPatch((x, 0.1), 0.22, 0.16, boxstyle="round,pad=.02", fc=color, ec="none", alpha=0.16))
    ax.text(x + 0.11, 0.18, text, ha="center", va="center", fontsize=6.4)


def arrow(ax: plt.Axes, x0: float, x1: float, y: float = 0.18) -> None:
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>", mutation_scale=9, lw=0.9, color=GRAY))


def grouped_lines(path: Path) -> list[tuple[str, np.ndarray, np.ndarray]]:
    table = rows(path)
    groups: dict[tuple[str, str], list[tuple[float, float]]] = {}
    for row in table:
        if row["kind"] != "line":
            continue
        groups.setdefault((row["series"], row["label"]), []).append((float(row["x"]), float(row["y"])))
    return [(f"curve {series}" if label.startswith("_") else label, np.asarray([p[0] for p in data]), np.asarray([p[1] for p in data])) for (series, label), data in groups.items()]


def fig1_visual_overview() -> None:
    fig, ax = plt.subplots(figsize=(7.15, 4.0))
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.94)
    image = plt.imread(FIG_DIR / "assets" / "fig1_afterglow_source.png")
    ax.imshow(image, extent=(0, 1, 0, 1), aspect="auto")
    ax.set_title("From a relativistic blast wave to an inspectable observer signal", fontsize=9)

    ax.annotate("central engine", (.126, .505), (.08, .63), arrowprops={"arrowstyle": "-", "color": GRAY}, ha="center")
    ax.annotate("relativistic ejecta", (.38, .50), (.34, .65), arrowprops={"arrowstyle": "-", "color": BLUE}, color=BLUE, ha="center")
    ax.annotate("reverse shock", (.56, .50), (.53, .73), arrowprops={"arrowstyle": "-", "color": RED}, color=RED, ha="center")
    ax.annotate("forward shock", (.66, .50), (.67, .78), arrowprops={"arrowstyle": "-", "color": TEAL}, color=TEAL, ha="center")

    inset = ax.inset_axes((.765, .285, .19, .43))
    inset.set_facecolor((1, 1, 1, .62))
    shell = [("ejecta", RED), ("RS", RED), ("contact", GRAY), ("FS", TEAL), ("medium", BLUE)]
    for index, (label, color) in enumerate(shell):
        x = .05 + .18 * index
        inset.add_patch(FancyBboxPatch((x, .16), .16, .62, boxstyle="round,pad=.01", fc=color, ec="white", alpha=.24))
        inset.text(x + .08, .48, label, rotation=90, ha="center", va="center", color=color, fontsize=6.2)
    inset.text(.5, .94, "magnified shell", ha="center", va="top", fontsize=6.5, fontweight="bold")
    inset.set(xlim=(0, 1), ylim=(0, 1))
    inset.axis("off")

    steps = [("1  store local state", .18, BLUE), ("2  transport particles + photons", .48, TEAL), ("3  project to the observer", .80, RED)]
    for text, x, color in steps:
        ax.text(x, .08, text, ha="center", va="center", color=color, fontsize=6.6, fontweight="bold",
                bbox={"boxstyle": "round,pad=.25", "fc": "white", "ec": color, "alpha": .86, "lw": .7})
    for x0, x1 in ((.29, .36), (.63, .69)):
        ax.add_patch(FancyArrowPatch((x0, .08), (x1, .08), arrowstyle="-|>", mutation_scale=8, color=GRAY, lw=.8))
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    save_pub(fig, "fig1_visual_overview")


def fig2_external_media() -> None:
    table = rows(DATA_DIR / "fig2_external_media_response.csv")
    media = (
        ("ism", "ISM", BLUE),
        ("wind_r2", r"wind $R^{-2}$", TEAL),
        ("explicit_jump", "Gaussian jump", RED),
        ("pion_lbv", r"PION $\eta$ Car CSM", GOLD),
    )
    fig, axes = plt.subplots(3, 4, figsize=(7.15, 5.0))
    fig.subplots_adjust(left=.07, right=.98, bottom=.10, top=.92, hspace=.34, wspace=.38)
    for column, (name, title, color) in enumerate(media):
        state = [row for row in table if row["medium"] == name and row["kind"] == "state"]
        flux = [row for row in table if row["medium"] == name and row["kind"] == "flux"]
        if not state or not flux:
            raise ValueError(f"fig2 source data missing state or flux rows for {name}")
        state_value = np.column_stack([values(state, key) for key in
                                       ("radius_cm", "density_cm3", "gamma_bulk", "bfield_g")])
        flux_value = np.column_stack([values(flux, key) for key in
                                      ("observer_time_s", "frequency_hz", "flux_cgs")])
        if not np.all(np.isfinite(state_value)) or not np.all(state_value > 0):
            raise ValueError(f"fig2 state contains non-finite or non-positive values for {name}")
        if not np.all(np.isfinite(flux_value)) or not np.all(flux_value > 0):
            raise ValueError(f"fig2 flux contains non-finite or non-positive values for {name}")
        radius = state_value[:, 0]
        if name == "pion_lbv":
            raw = [row for row in table if row["medium"] == name and row["kind"] == "profile_raw"]
            interface = [row for row in table if row["medium"] == name and row["kind"] == "profile_interface"]
            if len(raw) != 1968 or len(interface) != 96:
                raise ValueError("fig2 PION panel requires 1968 raw rows and 96 interface rows")
            raw_r, raw_n = values(raw, "radius_cm"), values(raw, "density_cm3")
            interface_r, interface_n = values(interface, "radius_cm"), values(interface, "density_cm3")
            axes[0, column].loglog(raw_r, raw_n, color=GRAY, lw=.65, alpha=.8, label="PION raw (2048)")
            axes[0, column].loglog(
                interface_r, interface_n, color=color, lw=1.35, ls="--", marker="o",
                markevery=8, ms=1.8, label="ASGARD interface (96)",
            )
            boundary_r = float(np.max(radius))
            for axis in axes[:2, column]:
                axis.axvspan(boundary_r, raw_r[-1], color=LIGHT, alpha=.65, zorder=0)
                axis.axvline(boundary_r, color=GRAY, ls=":", lw=.8)
                axis.set_xlim(raw_r[0], raw_r[-1])
            boundary_t = float(np.max(flux_value[:, 0]))
            axes[2, column].axvline(boundary_t, color=GRAY, ls=":", lw=.8)
            axes[2, column].axvspan(boundary_t, 2.5 * boundary_t, color=LIGHT, alpha=.65, zorder=0)
            axes[2, column].set_xlim(flux_value[0, 0], 2.5 * boundary_t)
            axes[0, column].text(
                .97, .96, "relativistic-solver\nboundary", transform=axes[0, column].transAxes,
                ha="right", va="top", color=GRAY, fontsize=4.8,
            )
            axes[0, column].legend(fontsize=4.6, loc="lower left")
        else:
            axes[0, column].loglog(radius, state_value[:, 1], color=color, lw=1.5)
        axes[1, column].loglog(radius, state_value[:, 2], color=BLUE, lw=1.3)
        field = axes[1, column].twinx()
        field.loglog(radius, state_value[:, 3], color=TEAL, ls="--", lw=1.1)
        field.set_ylabel(r"$B'$ (G)", color=TEAL, fontsize=6)
        field.tick_params(axis="y", labelcolor=TEAL, labelsize=5)
        axes[2, column].loglog(flux_value[:, 0], flux_value[:, 2], color=color, lw=1.4)
        panel(axes[0, column], chr(ord("a") + column))
        panel(axes[1, column], chr(ord("e") + column))
        panel(axes[2, column], chr(ord("i") + column))
        axes[0, column].set_title(title, fontsize=7.2)
        axes[0, column].set_ylabel(r"$n_{\rm p,eq}$ (cm$^{-3}$)" if column == 0 else "")
        axes[1, column].set_ylabel(r"$\Gamma$" if column == 0 else "")
        axes[2, column].set_ylabel(r"$F_\nu$" if column == 0 else "")
        axes[1, column].set_xlabel(r"$R$ (cm)")
        axes[2, column].set_xlabel(r"$t_{\rm obs}$ (s)")
        for axis in axes[:, column]:
            axis.grid(color=LIGHT, lw=.5, which="both")
        axes[0, column].text(.04, .08, r"$1\ n(R)$", transform=axes[0, column].transAxes,
                             color=color, fontsize=5.5)
        axes[1, column].text(.04, .08, r"$2\ \Gamma(R),B'(R)$",
                             transform=axes[1, column].transAxes, color=BLUE, fontsize=5.5)
        axes[2, column].text(.04, .08, r"$3\ F_\nu(t)$",
                             transform=axes[2, column].transAxes, color=color, fontsize=5.5)
    save_pub(fig, "fig2_external_media")


def fig3_radius_state() -> None:
    state=[r for r in rows(DATA_DIR/"fig2_forward_api.csv") if r["kind"]=="dynamics_state"]; local=[r for r in rows(BENCH/"cross_code"/"complex_state.csv") if r["code"]=="ASGARD" and r["value"]]; cross=rows(BENCH/"cross_code"/"fs.csv"); spectrum=[r for r in rows(DATA_DIR/"fig2_forward_api.csv") if r["kind"]=="spectrum_t1e+05s"]
    fig,axes=plt.subplots(2,3,figsize=(7.15,4.25)); fig.subplots_adjust(left=.08,right=.99,bottom=.12,top=.91,wspace=.42,hspace=.50)
    radius,gamma,field=values(state,"radius_cm"),values(state,"Gamma"),values(state,"B_comoving_G"); axes[0,0].loglog(radius,gamma,color=BLUE,lw=1.7); twin=axes[0,0].twinx(); twin.loglog(radius,field,color=TEAL,lw=1.4); twin.set_ylabel(r"$B'$ (G)",color=TEAL); panel(axes[0,0],"a"); finish(axes[0,0],r"$R$ (cm)",r"$\Gamma$","Radius-ordered state")
    for metric,color,label in (("electron_number",BLUE,r"$N_e$"),("electron_energy_erg",VIOLET,r"$E_e$")):
        use=sorted([r for r in local if r["value_name"]==metric],key=lambda r:float(r["radius_cm"])); axes[0,1].loglog(values(use,"radius_cm"),values(use,"value"),color=color,lw=1.4,label=label)
    panel(axes[0,1],"b"); finish(axes[0,1],r"$R$ (cm)","stored state","Particles and energy"); axes[0,1].legend(fontsize=5.5)
    chi_value=min(float(r["chi"]) for r in local if r["metric_group"]=="chi_projection" and r["chi"])
    for metric,color,line,label in (("max_tau_syn_chi",GRAY,"-",rf"$\max\tau_\nu$, $\chi={chi_value:.2f}$"),("max_seed_syn_chi",TEAL,"--",rf"seed source, $\chi={chi_value:.2f}$")):
        use=sorted([r for r in local if r["value_name"]==metric and r["chi"] and float(r["chi"])==chi_value],key=lambda r:float(r["radius_cm"])); axes[0,2].loglog(values(use,"radius_cm"),values(use,"value"),color=color,ls=line,lw=1.3,label=label)
    panel(axes[0,2],"c"); finish(axes[0,2],r"$R$ (cm)","stored transfer state","One downstream trajectory"); axes[0,2].legend(fontsize=5)
    styles=(("ASGARD",BLUE,"-","o"),("afterglowpy",TEAL,"--","s"),("jetsimpy",VIOLET,"-.","^"),("VegasAfterglow",RED,":","D"),("PyBlastAfterglowMag",GOLD,(0,(4,1.3)),"v")); ratios=[]
    for code,color,line,marker in styles:
        sub=[r for r in cross if r["code"]==code and r["value_name"]=="flux_density_cgs" and float(r["frequency_hz"])==1e14]; axes[1,0].loglog(values(sub,"observer_time"),values(sub,"value"),color=color,ls=line,marker=marker,markevery=4,ms=2,lw=1.1,label=code)
        if code!="ASGARD": ratios.append((code,color,line,marker,[r for r in sub if r["asgard_ratio"]]))
    left=max(min(values(group,"observer_time")) for _,_,_,_,group in ratios); right=min(max(values(group,"observer_time")) for _,_,_,_,group in ratios); axes[1,0].axvspan(left,right,color=BLUE,alpha=.055,zorder=0); panel(axes[1,0],"d"); finish(axes[1,0],r"$t_{\rm obs}$ (s)",r"$F_\nu$","Matched five-code FS light curve"); axes[1,0].legend(fontsize=4.5,ncol=2)
    axes[1,1].axvspan(left,right,color=BLUE,alpha=.08,zorder=0)
    for code,color,line,marker,group in ratios: axes[1,1].loglog(values(group,"observer_time"),values(group,"asgard_ratio"),color=color,ls=line,marker=marker,markevery=4,ms=2,lw=1.1)
    axes[1,1].axhline(1,color=GRAY,lw=.8); axes[1,1].text(.03,.06,"shared finite positive domain",transform=axes[1,1].transAxes,color=BLUE,fontsize=5.2); panel(axes[1,1],"e"); finish(axes[1,1],r"$t_{\rm obs}$ (s)","code / ASGARD","Common-domain ratio")
    axes[1,2].loglog(values(spectrum,"frequency_hz"),values(spectrum,"frequency_hz")*values(spectrum,"total_fnu_cgs"),color=GRAY,lw=1.5,label="total"); axes[1,2].loglog(values(spectrum,"frequency_hz"),values(spectrum,"frequency_hz")*values(spectrum,"fs_synch_fnu_cgs"),color=TEAL,ls="--",lw=1.2,label="FS synchrotron"); axes[1,2].loglog(values(spectrum,"frequency_hz"),values(spectrum,"frequency_hz")*values(spectrum,"fs_ssc_fnu_cgs"),color=VIOLET,ls="-.",lw=1.2,label="FS SSC"); panel(axes[1,2],"f"); finish(axes[1,2],r"$\nu$ (Hz)",r"$\nu F_\nu$","Observer SED at $10^5$ s"); axes[1,2].legend(fontsize=5.2)
    save_pub(fig,"fig3_radius_state")


def fig4_jet_profiles() -> None:
    theta=np.linspace(0,.4,300); core=.08; power=np.ones_like(theta); wing=theta>core; power[wing]=(theta[wing]/core)**-2
    profiles=(("top-hat",(theta<=core).astype(float),BLUE,"-"),("Gaussian",np.exp(-.5*(theta/core)**2),TEAL,"--"),("power law",power,VIOLET,"-."))
    viewing=[row for row in rows(BENCH/"angular_reference"/"fullref300x48_vs_dominant-region-ioka-time-v1_20x15_g120_tc0p08_lightcurves.csv") if float(row["frequency_hz"])==1e14]
    cross=[row for row in rows(BENCH/"cross_code"/"rs_ssc_geometry.csv") if row["scenario"]=="gaussian_off_axis_fs_synch" and row["value_name"]=="fs_synch_flux_density_cgs" and float(row["frequency_hz"])==5e14]
    fig,axes=plt.subplots(2,3,figsize=(7.15,4.2)); fig.subplots_adjust(left=.07,right=.99,bottom=.13,top=.91,wspace=.42,hspace=.48)
    for name,profile,color,line in profiles:
        axes[0,0].plot(theta/core,profile,color=color,ls=line,lw=1.5,label=name); axes[0,1].semilogy(theta/core,1+299*profile,color=color,ls=line,lw=1.5,label=name)
    panel(axes[0,0],"a"); finish(axes[0,0],r"$\theta/\theta_c$",r"$(dE/d\Omega)/(dE/d\Omega)_0$","Angular energy"); axes[0,0].legend(fontsize=5)
    panel(axes[0,1],"b"); finish(axes[0,1],r"$\theta/\theta_c$",r"$\Gamma(\theta)$","Lorentz factor")
    axes[0,2].axis("off"); axes[0,2].legend(handles=[Line2D([],[],color=color,ls=line,lw=1.5,label=name) for name,_,color,line in profiles],loc="center",fontsize=6)
    for index,ratio in enumerate(sorted({float(row["theta_obs_over_theta_c"]) for row in viewing})):
        use=[row for row in viewing if float(row["theta_obs_over_theta_c"])==ratio]; axes[1,0].loglog(values(use,"time_s"),values(use,"reference"),color=COLORS[index],ls=("-","--","-.",":",(0,(4,1.3)))[index],lw=1.1,label=f"{ratio:g} theta_c")
    panel(axes[1,0],"c"); finish(axes[1,0],r"$t_{\rm obs}$ (s)",r"$F_\nu$","Viewing-angle sequence"); axes[1,0].legend(fontsize=4.7,ncol=2)
    styles=(("ASGARD",BLUE,"-"),("afterglowpy",TEAL,"--"),("jetsimpy",VIOLET,"-."),("VegasAfterglow",RED,":"),("PyBlastAfterglowMag",GOLD,(0,(4,1.3))))
    for code,color,line in styles:
        use=[row for row in cross if row["code"]==code]; axes[1,1].loglog(values(use,"observer_time"),values(use,"value"),color=color,ls=line,lw=1.1,label=code)
        if code!="ASGARD": axes[1,2].loglog(values(use,"observer_time"),values(use,"asgard_ratio"),color=color,ls=line,lw=1.1)
    panel(axes[1,1],"d"); finish(axes[1,1],r"$t_{\rm obs}$ (s)",r"$F_\nu$","Five-code Gaussian"); axes[1,1].legend(fontsize=4.2,ncol=2)
    axes[1,2].axhline(1,color=BLUE,lw=.8); panel(axes[1,2],"e"); finish(axes[1,2],r"$t_{\rm obs}$ (s)","code / ASGARD","Matched ratio")
    save_pub(fig,"fig4_jet_profiles")


def fig5_electron_transport() -> None:
    table=rows(DATA_DIR/"fig3_electron_transport.csv")
    spec=[r for r in table if r["kind"]=="electron_spectrum"]; budget=[r for r in table if r["kind"]=="shell_budget"]
    radius=np.unique(values(spec,"radius_cm")); gamma=np.unique(values(spec,"gamma_e")); electron=np.asarray([float(r["dN_dgamma_e"]) for r in spec]).reshape(gamma.size,radius.size)
    fig,axes=plt.subplots(1,3,figsize=(7.15,2.7)); fig.subplots_adjust(left=.07,right=.99,bottom=.22,top=.84,wspace=.45)
    image=axes[0].pcolormesh(np.log10(radius),np.log10(gamma),np.ma.log10(np.ma.masked_less_equal(electron,0.0)),shading="auto",cmap="viridis")
    panel(axes[0],"a"); axes[0].set(xlabel=r"$\log R$ (cm)",ylabel=r"$\log\gamma_e$",title="Electron transport history"); fig.colorbar(image,ax=axes[0],fraction=.045,pad=.02,label=r"$\log(dN_e/d\gamma_e)$")
    for idx,color in zip(np.linspace(0,radius.size-1,4,dtype=int),COLORS): axes[1].loglog(gamma,gamma**2*electron[:,idx],color=color,lw=1.2,label=rf"$R={radius[idx]:.1e}$")
    panel(axes[1],"b"); finish(axes[1],r"$\gamma_e$",r"$\gamma_e^2dN_e/d\gamma_e$","Injection, cooling, aging"); axes[1].legend(fontsize=5)
    rb=values(budget,"radius_cm"); axes[2].loglog(rb,values(budget,"electron_number"),color=BLUE,lw=1.5,label=r"$N_e$"); twin=axes[2].twinx(); twin.loglog(rb,values(budget,"electron_energy_erg"),color=VIOLET,lw=1.3); twin.set_ylabel(r"$E_e$ (erg)",color=VIOLET)
    panel(axes[2],"c"); finish(axes[2],r"$R$ (cm)",r"$N_e$","Number and energy budgets")
    save_pub(fig,"fig5_electron_transport")


def fig6_electron_solvers() -> None:
    data = BENCH / "vegas_comparison" / "baseline"
    spectra = rows(data / "electron_solver_spectra.csv")
    observables = rows(data / "electron_solver_observables.csv")
    solvers = ("fullhide_1d", "fullhide_2d", "dg_1d", "slc1_1d", "charint_1d", "charint_2d")
    style = {
        "fullhide_1d": (BLUE, 1.7, 1.0, "-"),
        "fullhide_2d": (TEAL, 1.5, 1.0, ":"),
        "dg_1d": (VIOLET, 1.7, 1.0, "--"),
        "slc1_1d": (RED, .85, .55, "-."),
        "charint_1d": (GOLD, .85, .55, (0, (4, 1.4))),
        "charint_2d": (GRAY, .85, .55, (0, (1, 1.2))),
    }
    epoch = 1.0e5
    fig, axes = plt.subplots(2, 3, figsize=(7.15, 4.5))
    fig.subplots_adjust(left=.08, right=.99, bottom=.12, top=.91, wspace=.45, hspace=.52)

    for solver in solvers:
        color, width, alpha, line = style[solver]
        spec = [row for row in spectra if row["kind"] == "spectrum" and row["solver"] == solver and float(row["requested_time_s"]) == epoch]
        axes[0, 0].loglog(values(spec, "gamma_e"), values(spec, "value"), color=color, lw=width, alpha=alpha, ls=line, label=solver)
        budget = [row for row in spectra if row["kind"] == "budget" and row["solver"] == solver]
        axes[0, 1].loglog(values(budget, "requested_time_s"), values(budget, "electron_number"), color=color, lw=width, alpha=alpha, ls=line)
        axes[0, 2].loglog(values(budget, "requested_time_s"), values(budget, "electron_energy_erg"), color=color, lw=width, alpha=alpha, ls=line)
        sed = [row for row in observables if row["kind"] == "sed" and row["solver"] == solver and float(row["time_s"]) == epoch]
        axes[1, 0].loglog(values(sed, "frequency_hz"), values(sed, "value"), color=color, lw=width, alpha=alpha, ls=line)
        lightcurve = [row for row in observables if row["kind"] == "lightcurve" and row["solver"] == solver]
        axes[1, 1].loglog(values(lightcurve, "time_s"), values(lightcurve, "value"), color=color, lw=width, alpha=alpha, ls=line)

    vegas_sed = [row for row in observables if row["kind"] == "sed" and row["solver"] == "VegasAfterglow" and float(row["time_s"]) == epoch]
    vegas_lightcurve = [row for row in observables if row["kind"] == "lightcurve" and row["solver"] == "VegasAfterglow"]
    axes[1, 0].loglog(values(vegas_sed, "frequency_hz"), values(vegas_sed, "value"), color=GRAY, lw=1.0, ls="--")
    axes[1, 1].loglog(values(vegas_lightcurve, "time_s"), values(vegas_lightcurve, "value"), color=GRAY, lw=1.0, ls="--")

    reference=[row for row in spectra if row["kind"]=="spectrum" and row["solver"]=="fullhide_1d" and float(row["requested_time_s"])==epoch]
    ref_gamma=values(reference,"gamma_e"); ref_value=values(reference,"value"); ref_positive=ref_value>0; ref_gamma=ref_gamma[ref_positive]; ref_value=ref_value[ref_positive]
    for solver in solvers[1:]:
        color,width,alpha,line=style[solver]
        spec=[row for row in spectra if row["kind"]=="spectrum" and row["solver"]==solver and float(row["requested_time_s"])==epoch]
        solver_gamma=values(spec,"gamma_e"); solver_value=values(spec,"value")
        common=(solver_gamma>=ref_gamma.min())&(solver_gamma<=ref_gamma.max())&(solver_value>0)
        ref_on=np.exp(np.interp(np.log(solver_gamma[common]),np.log(ref_gamma),np.log(ref_value)))
        axes[1,2].semilogx(solver_gamma[common],np.log10(solver_value[common]/ref_on),color=color,lw=width,alpha=alpha,ls=line)
    axes[1,2].axhline(0.0,color=BLUE,lw=.8)

    labels = (
        (axes[0, 0], "a", r"$\gamma_e$", r"$dN_e/d\gamma_e$", r"Electron spectra at $10^5$ s"),
        (axes[0, 1], "b", r"$t_{\rm obs}$ (s)", r"$N_e$", "Electron-number budget"),
        (axes[0, 2], "c", r"$t_{\rm obs}$ (s)", r"$E_e$ (erg)", "Electron-energy budget"),
        (axes[1, 0], "d", r"$\nu$ (Hz)", r"$\nu F_\nu$", r"SED at $10^5$ s"),
        (axes[1, 1], "e", r"$t_{\rm obs}$ (s)", r"$F_\nu$", r"Light curve at $10^{14}$ Hz"),
        (axes[1, 2], "f", r"$\gamma_e$", r"$\log_{10}(N_e/N_{e,\rm fullhide})$", "Electron-spectrum difference"),
    )
    for axis, letter, xlabel, ylabel, title in labels:
        panel(axis, letter)
        finish(axis, xlabel, ylabel, title)
    axes[0, 0].legend(fontsize=4.8, ncol=2)
    save_pub(fig,"fig6_electron_solvers")


def fig7_radiation() -> None:
    table=[r for r in rows(DATA_DIR/"fig2_forward_api.csv") if r["kind"]=="spectrum_t1e+05s"]
    nu=values(table,"frequency_hz"); syn=values(table,"fs_synch_fnu_cgs"); ssc=values(table,"fs_ssc_fnu_cgs"); total=values(table,"total_fnu_cgs")
    fig,axes=plt.subplots(2,2,figsize=(7.15,4.15)); fig.subplots_adjust(left=.08,right=.98,bottom=.13,top=.91,wspace=.40,hspace=.50)
    axes[0,0].loglog(nu,syn,color=TEAL,lw=1.7,label="synchrotron"); axes[0,0].loglog(nu,ssc,color=VIOLET,lw=1.7,label="SSC"); axes[0,0].loglog(nu,total,color=GRAY,lw=.9,label="total")
    panel(axes[0,0],"a"); finish(axes[0,0],r"$\nu$ (Hz)",r"$F_\nu$","Electron spectrum to radiation"); axes[0,0].legend()
    frac=np.divide(ssc,total,out=np.zeros_like(ssc),where=total>0); axes[0,1].semilogx(nu,frac,color=VIOLET,lw=1.6); axes[0,1].fill_between(nu,0,frac,color=VIOLET,alpha=.12)
    panel(axes[0,1],"b"); finish(axes[0,1],r"$\nu$ (Hz)",r"$F_\nu^{SSC}/F_\nu$","Component-dominance map"); axes[0,1].set_ylim(0,1)

    state=rows(BENCH/"cross_code"/"complex_state.csv")
    tau_rows=[row for row in state if row["metric_group"]=="chi_projection" and row["value_name"]=="max_tau_syn_chi" and row["value"]]
    chi=sorted({float(row["chi"]) for row in tau_rows})
    groups=[sorted([row for row in tau_rows if float(row["chi"])==number],key=lambda row:float(row["radius_cm"])) for number in chi]
    radius_grid=np.asarray([[float(row["radius_cm"]) for row in group] for group in groups])
    matrix=np.asarray([[float(row["value"]) for row in group] for group in groups])
    chi_grid=np.broadcast_to(np.asarray(chi)[:,None],matrix.shape)
    image=axes[1,0].pcolormesh(np.log10(radius_grid),chi_grid,np.ma.log10(np.ma.masked_less_equal(matrix,0.0)),shading="nearest",cmap="magma")
    panel(axes[1,0],"c"); axes[1,0].set(xlabel=r"$\log_{10}R$ (cm)",ylabel=r"$\chi$",title=r"Stored $\chi$-resolved SSA depth"); fig.colorbar(image,ax=axes[1,0],fraction=.045,pad=.02,label=r"$\log_{10}\max\tau_\nu$")

    transfer=[row for row in rows(DATA_DIR/"figA1_hadronic_thresholds.csv") if row["kind"]=="survival"]
    optical_depth=values(transfer,"target_photon_ev")
    axes[1,1].loglog(optical_depth,values(transfer,"pgamma_gamma_p_threshold"),color=GRAY,lw=1.6,label=r"$e^{-\tau}$")
    axes[1,1].loglog(optical_depth,values(transfer,"bh_gamma_p_threshold"),color=RED,lw=1.6,label=r"$(1-e^{-\tau})/\tau$")
    panel(axes[1,1],"d"); finish(axes[1,1],r"shell optical depth $\tau$","survival / cell transfer",r"$\gamma\gamma$ sink semantics"); axes[1,1].legend()
    save_pub(fig,"fig7_radiation")


def fig8_primary_reverse_shock() -> None:
    comp=rows(BENCH/"cross_code"/"rs_ssc_geometry.csv")
    fig,axes=plt.subplots(1,3,figsize=(7.15,2.75)); fig.subplots_adjust(left=.06,right=.99,bottom=.22,top=.84,wspace=.44)
    shell=(("ejecta",RED),("RS",RED),("contact",GRAY),("FS",TEAL),("medium",BLUE))
    for index,(label,color) in enumerate(shell):
        x=.04+.18*index; axes[0].add_patch(FancyBboxPatch((x,.24),.16,.48,boxstyle="round,pad=.01",fc=color,ec="white",alpha=.24)); axes[0].text(x+.08,.48,label,rotation=90,ha="center",va="center",color=color)
    axes[0].set(xlim=(0,1),ylim=(0,1)); axes[0].axis("off"); panel(axes[0],"a"); axes[0].set_title("Primary ejecta shock structure")

    component_style=(("fs_sync_flux_density_cgs",TEAL,"FS synchrotron"),("fs_ssc_flux_density_cgs",VIOLET,"FS SSC"),("rs_sync_flux_density_cgs",RED,"RS synchrotron"),("rs_ssc_flux_density_cgs",GOLD,"RS SSC"))
    for name,color,label in component_style:
        use=sorted([row for row in comp if row["scenario"]=="rs_ssc_matched_tophat" and row["code"]=="ASGARD" and row["value_name"]==name and row["value"] and float(row["frequency_hz"])==5.0e14],key=lambda row:float(row["observer_time"]))
        axes[1].loglog(values(use,"observer_time"),values(use,"value"),color=color,lw=1.3,label=label)
    panel(axes[1],"b"); finish(axes[1],r"$t_{\rm obs}$ (s)",r"$F_\nu$","FS--RS component hand-off"); axes[1].legend(fontsize=4.8)

    for name,color,line,label in (("rs_sync_flux_density_cgs",RED,"-","Vegas RS synchrotron"),("rs_ssc_flux_density_cgs",GOLD,"--","Vegas RS SSC")):
        use=sorted([r for r in comp if r["scenario"]=="rs_ssc_matched_tophat" and r["code"]=="VegasAfterglow" and r["value_name"]==name and r["asgard_ratio"] and float(r["frequency_hz"])==5.0e14],key=lambda r:float(r["observer_time"]))
        axes[2].loglog(values(use,"observer_time"),values(use,"asgard_ratio"),color=color,ls=line,lw=1.4,label=label)
    axes[2].axhline(1,color=GRAY,lw=.8); panel(axes[2],"c"); finish(axes[2],r"$t_{\rm obs}$ (s)","Vegas / ASGARD","Matched RS component ratios"); axes[2].legend(fontsize=4.8)
    save_pub(fig,"fig8_primary_reverse_shock")


def fig9_magnetized_reverse_shock() -> None:
    fig,axes=plt.subplots(2,3,figsize=(7.15,4.2)); fig.subplots_adjust(left=.08,right=.99,bottom=.12,top=.92,wspace=.45,hspace=.48)
    for medium,ls in (("ism","-"),("wind","--")):
        summary=rows(BENCH/"magnetized_sigma"/medium/"sigma_scan_summary.csv"); sig=values(summary,"sigma"); xp=sig
        axes[0,0].plot(xp,values(summary,"max_B3_G"),ls=ls,color=RED,lw=1.4,label=medium)
        axes[0,1].plot(xp,values(summary,"max_gamma34"),ls=ls,color=BLUE,lw=1.4)
        lc=[r for r in rows(BENCH/"magnetized_sigma"/medium/"sigma_scan_lightcurve_summary.csv") if float(r["nu_hz"])==1e14]
        axes[0,2].plot(xp,values(lc,"rs_to_total_at_total_peak"),ls=ls,color=VIOLET,lw=1.4)
        axes[1,0].plot(xp,values(summary,"t_cross_s"),ls=ls,color=GOLD,lw=1.4)
        raw=rows(BENCH/"magnetized_sigma"/medium/"sigma_scan_sed_raw.csv")
        for sigma,color in zip((0.0,.1,1.0,100.0),COLORS):
            sed=sorted([row for row in raw if float(row["sigma"])==sigma and float(row["time_s"])==1.0e3],key=lambda row:float(row["frequency_hz"]))
            axes[1,1 if medium=="ism" else 2].loglog(values(sed,"frequency_hz"),values(sed,"frequency_hz")*values(sed,"total_fnu_cgs"),color=color,lw=1.1,label=rf"$\sigma={sigma:g}$")
    for axis in (axes[0,0],axes[0,1],axes[0,2],axes[1,0]): axis.set_xscale("symlog",linthresh=1.0e-4)
    axes[0,0].set_yscale("log"); axes[1,0].set_yscale("log")
    labels=(
        (axes[0,0],"a",r"$\sigma$",r"max $B'_3$ (G)","Compressed magnetic field"),
        (axes[0,1],"b",r"$\sigma$",r"max $\gamma_{34}$","Shock strength"),
        (axes[0,2],"c",r"$\sigma$","RS fraction","Optical peak contribution"),
        (axes[1,0],"d",r"$\sigma$",r"$t_\times$ (s)","Crossing time"),
        (axes[1,1],"e",r"$\nu$ (Hz)",r"$\nu F_\nu$","ISM SED at $10^3$ s"),
        (axes[1,2],"f",r"$\nu$ (Hz)",r"$\nu F_\nu$","Wind SED at $10^3$ s"),
    )
    for ax,label,xlabel,ylabel,title in labels: panel(ax,label); finish(ax,xlabel,ylabel,title)
    axes[0,0].legend()
    axes[1,1].legend(fontsize=4.8,ncol=2)
    save_pub(fig,"fig9_magnetized_reverse_shock")


def density_structure_figure(stem: str, title: str) -> None:
    data=BENCH/"density_structure"; profile=rows(data/f"{stem}_profile.csv"); dynamics=rows(data/f"{stem}_dynamics.csv"); events=rows(data/f"{stem}_events.csv"); energy=rows(data/f"{stem}_secondary_rs_energy.csv"); flux=rows(data/f"{stem}_flux.csv")
    fig,axes=plt.subplots(2,2,figsize=(7.15,4.0)); fig.subplots_adjust(left=.08,right=.99,bottom=.13,top=.91,wspace=.40,hspace=.48)
    axes[0,0].loglog(values(profile,"radius_cm"),values(profile,"density_cm3"),color=BLUE,lw=1.5)
    for event in events: axes[0,0].axvline(float(event["input_center_cm"]),color=RED if event["event_active"]=="True" else GRAY,ls=":" if event["event_active"]=="True" else "--",lw=.7)
    panel(axes[0,0],"a"); finish(axes[0,0],r"$R$ (cm)",r"$n(R)$ (cm$^{-3}$)",title)
    for index,event in enumerate(events):
        branch=[row for row in dynamics if int(row["event_index"])==index and float(row["branch_internal_energy_erg"])>0]
        if branch: axes[0,1].loglog(values(branch,"radius_cm"),values(branch,"branch_internal_energy_erg"),color=COLORS[index],ls=("-","--","-.")[index],lw=1.2,label=f"event {index}")
        if event["event_active"]=="True": axes[0,1].axvspan(float(event["start_radius_cm"]),float(event["shock_end_radius_cm"]),color=COLORS[index],alpha=.08)
        else: axes[0,1].text(.55,.12+index*.1,f"event {index}: no shock",transform=axes[0,1].transAxes,color=GRAY,fontsize=5.5)
    panel(axes[0,1],"b"); finish(axes[0,1],r"$R$ (cm)","branch internal energy (erg)","Tracked event state"); axes[0,1].legend(fontsize=5)
    light=[row for row in flux if float(row["band_hz"])==1e14]
    for key,color,line,label in (("forward_sync_flux_cgs",TEAL,"-","FS"),("reverse_sync_flux_cgs",RED,"--","secondary RS"),("total_flux_cgs",GRAY,"-.","total")): axes[1,0].loglog(values(light,"time_s"),values(light,key),color=color,ls=line,lw=1.2,label=label)
    panel(axes[1,0],"c"); finish(axes[1,0],r"$t_{\rm obs}$ (s)",r"$F_\nu$",r"Components at $10^{14}$ Hz"); axes[1,0].legend(fontsize=5)
    event_energy=[row for row in energy if row["event_index"]!="total"]; x=np.arange(len(event_energy)); dissipated=values(event_energy,"secondary_rs_dissipated_energy_erg"); injected=values(event_energy,"secondary_rs_electron_injected_energy_erg")
    axes[1,1].bar(x-.16,dissipated,.32,color=RED,alpha=.8,label="dissipated"); axes[1,1].bar(x+.16,injected,.32,color=GOLD,alpha=.8,label="electron injection")
    positive=np.concatenate((dissipated[dissipated>0],injected[injected>0])); axes[1,1].set_yscale("log"); axes[1,1].set_ylim(positive.min()/3,positive.max()*3)
    for index,row in enumerate(event_energy):
        if float(row["secondary_rs_dissipated_energy_erg"])==0: axes[1,1].text(index,positive.min()/2,"no event",ha="center",va="top",fontsize=5,color=GRAY)
    axes[1,1].set_xticks(x,[f"event {row['event_index']}" for row in event_energy]); panel(axes[1,1],"d"); finish(axes[1,1],"","energy (erg)","Secondary-shock energy"); axes[1,1].legend(fontsize=5)
    save_pub(fig,stem)


def fig10_wind_termination_shell() -> None: density_structure_figure("fig10_wind_termination_shell","Wind termination + finite shell")
def fig11_density_jumps() -> None:
    data=BENCH/"density_jump"; meta=json.loads((data/"metadata.json").read_text(encoding="utf-8"))
    events=rows(data/"triple_density_jump_secondary_rs_events.csv"); energy=rows(data/"triple_density_jump_secondary_rs_energy.csv"); flux=rows(data/"triple_density_jump_flux.csv")
    centers=np.asarray(meta["jump_radii_cm"]); factor=float(meta["jump_factor"]); width=float(meta["jump_width_relative"])
    radius=np.logspace(np.log10(centers.min())-1,np.log10(centers.max())+1,480); density=np.ones_like(radius)
    for center in centers: density+=(factor-1)*np.exp(-.5*((radius-center)/(width*center))**2)
    fig,axes=plt.subplots(2,2,figsize=(7.15,4.0)); fig.subplots_adjust(left=.08,right=.99,bottom=.13,top=.91,wspace=.38,hspace=.48)
    axes[0,0].loglog(radius,density,color=BLUE,lw=1.6)
    for center in centers: axes[0,0].axvline(center,color=BLUE,lw=.55,ls=":",alpha=.6)
    panel(axes[0,0],"a"); finish(axes[0,0],r"$R$ (cm)",r"$n(R)/n_0$","Tracked triple-jump input")
    active=[row for row in events if row["event_active"]=="True"]; left=min(float(row["start_tobs_axis_s"]) for row in active); right=max(float(row["shock_end_tobs_axis_s"]) for row in active)
    for index,event in enumerate(events):
        if event["event_active"]=="True":
            start=float(event["start_tobs_axis_s"]); end=float(event["shock_end_tobs_axis_s"]); axes[0,1].barh(index,end-start,left=start,color=RED,alpha=.7)
        else: axes[0,1].text(np.sqrt(left*right),index,"no event",ha="center",va="center",color=GRAY,fontsize=6)
    axes[0,1].set_xscale("log"); axes[0,1].set_yticks(range(len(events)),[f"event {index}" for index in range(len(events))]); panel(axes[0,1],"b"); finish(axes[0,1],r"$t_{\rm obs}$ (s)","","Secondary-RS windows")
    deposit=[row for row in energy if row["event_index"]!="total"]; x=np.arange(len(deposit)); dissipated=values(deposit,"secondary_rs_dissipated_energy_erg"); injected=values(deposit,"secondary_rs_electron_injected_energy_erg")
    axes[1,0].bar(x-.16,dissipated,.32,color=RED,label="dissipated"); axes[1,0].bar(x+.16,injected,.32,color=GOLD,label="electron injection")
    positive=np.concatenate((dissipated[dissipated>0],injected[injected>0])); axes[1,0].set_yscale("log"); axes[1,0].set_ylim(positive.min()/3,positive.max()*3); axes[1,0].set_xticks(x,[f"event {row['event_index']}" for row in deposit])
    for index,row in enumerate(deposit):
        if float(row["secondary_rs_dissipated_energy_erg"])==0: axes[1,0].text(index,positive.min()/2,"no event",ha="center",va="top",fontsize=5,color=GRAY)
    panel(axes[1,0],"c"); finish(axes[1,0],"","energy (erg)","Secondary energy deposition"); axes[1,0].legend(fontsize=5.2)
    light=[row for row in flux if float(row["band_hz"])==1.0e14]
    for key,color,line,label in (("jump_forward_sync",TEAL,"-","FS"),("jump_reverse_sync",RED,"--","secondary RS"),("jump_total",GRAY,"-.","total")): axes[1,1].loglog(values(light,"time_s"),values(light,key),color=color,ls=line,lw=1.2,label=label)
    panel(axes[1,1],"d"); finish(axes[1,1],r"$t_{\rm obs}$ (s)",r"$F_\nu$",r"Response at $10^{14}$ Hz"); axes[1,1].legend(fontsize=5.2)
    save_pub(fig,"fig11_density_jumps")
def fig12_smooth_clumpy_medium() -> None: density_structure_figure("fig12_smooth_clumpy_medium","Smooth tabulated clumps")


def fig13_hadronic_pair_state() -> None:
    state = [row for row in rows(BENCH / "cross_code" / "complex_state.csv") if row["code"] == "ASGARD" and row["value"]]
    threshold_table = rows(DATA_DIR / "figA1_hadronic_thresholds.csv")
    threshold = [row for row in threshold_table if row["kind"] == "threshold"]
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 4.15))
    fig.subplots_adjust(left=.08, right=.98, bottom=.13, top=.91, wspace=.42, hspace=.50)

    lines=[]
    for metric,color,label in (("proton_energy_erg",GOLD,r"legacy $E_p$"),("bh_pair_energy_erg","#D9A441",r"AM3 $E_{e^\pm}^{\rm BH}$")):
        use=sorted([row for row in state if row["value_name"]==metric],key=lambda row:float(row["radius_cm"]))
        lines.extend(axes[0,0].loglog(values(use,"radius_cm"),values(use,"value"),color=color,lw=1.5,label=label))
    number_axis=axes[0,0].twinx()
    for metric,color,label in (("proton_number",GOLD,r"legacy $N_p$"),("bh_pair_number","#D9A441",r"AM3 $N_{e^\pm}^{\rm BH}$")):
        use=sorted([row for row in state if row["value_name"]==metric],key=lambda row:float(row["radius_cm"]))
        lines.extend(number_axis.loglog(values(use,"radius_cm"),values(use,"value"),color=color,ls="--",lw=1.1,label=label))
    number_axis.set_ylabel("stored particle number")
    panel(axes[0,0],"a"); finish(axes[0,0],r"$R$ (cm)","stored energy (erg)","Tracked proton and pair states")
    axes[0,0].legend(lines,[line.get_label() for line in lines],fontsize=4.7,ncol=2)

    for metric,color,line,label in (("hadronic_synch_peak_cgs",GOLD,"-","proton synchrotron"),("pgamma_peak_cgs","#B07D00","--",r"$p\gamma$ peak"),("bethe_heitler_peak_cgs","#D9A441","-.","BH peak")):
        use = sorted([row for row in state if row["value_name"] == metric], key=lambda row: float(row["radius_cm"]))
        axes[0, 1].loglog(values(use, "radius_cm"), values(use, "value"), color=color, ls=line, lw=1.4, label=label)
    panel(axes[0, 1], "b")
    finish(axes[0, 1], r"$R$ (cm)", "stored peak (cgs)", "Activated radiative products")
    axes[0, 1].legend(fontsize=5.5)

    for metric,color,line,label in (("max_tau_pg",GOLD,"-",r"max $\tau_{p\gamma}$"),("max_tau_bh","#D9A441","--",r"max $\tau_{\rm BH}$")):
        use = sorted([row for row in state if row["value_name"] == metric], key=lambda row: float(row["radius_cm"]))
        axes[1, 0].loglog(values(use, "radius_cm"), values(use, "value"), color=color, ls=line, lw=1.4, label=label)
    survival = sorted([row for row in state if row["value_name"] == "min_pg_photon_survival"], key=lambda row: float(row["radius_cm"]))
    survival_axis = axes[1, 0].twinx()
    survival_axis.semilogx(values(survival, "radius_cm"), values(survival, "value"), color=GRAY, lw=1.3)
    survival_axis.set_ylabel("minimum photon survival", color=GRAY)
    panel(axes[1, 0], "c")
    finish(axes[1, 0], r"$R$ (cm)", "maximum optical depth", "Interaction depth and photon survival")
    axes[1, 0].legend(fontsize=5.5)

    eps = values(threshold, "target_photon_ev")
    axes[1, 1].loglog(eps, values(threshold, "pgamma_gamma_p_threshold"), color=GOLD, lw=1.5, label=r"$p\gamma$")
    axes[1, 1].loglog(eps, values(threshold, "bh_gamma_p_threshold"), color="#D9A441", ls="--", lw=1.5, label="BH")
    panel(axes[1, 1], "d")
    finish(axes[1, 1], r"target photon energy $\epsilon'$ (eV)", r"threshold $\gamma_p$", "Implemented interaction thresholds")
    axes[1, 1].legend()
    save_pub(fig,"fig13_hadronic_pair_state")


def fig14_eats_projection() -> None:
    conv=rows(BENCH/"chi_eats"/"chi_eats_convergence.csv")
    fig,axes=plt.subplots(1,3,figsize=(7.15,2.8),gridspec_kw={"width_ratios":[1.05,1,1]}); fig.subplots_adjust(left=.05,right=.99,bottom=.22,top=.84,wspace=.42)
    ax=axes[0]; ax.add_patch(Circle((.18,.48),.035,fc=GOLD)); ax.add_patch(Wedge((.18,.48),.65,-30,30,width=.012,fc=TEAL)); ax.plot([.18,.85],[.48,.48],color=GRAY,lw=.8); ax.plot([.18,.72],[.48,.75],color=BLUE,lw=1.2); ax.scatter([.88],[.48],marker=(3,0,-90),s=70,color=GRAY); ax.text(.45,.63,r"$R,\mu,\chi$",color=BLUE); ax.text(.72,.42,r"$\theta_{\rm obs}$"); ax.text(.5,.12,r"$(1+z)[t_{\rm on}(R)+R(1-\mu)/c]$",ha="center",fontsize=6.8); ax.axis("off"); panel(ax,"a"); ax.set_title("Equal-arrival-time geometry")
    for scan,color,marker,line in (("num_chi",BLUE,"o","-"),("num_theta",TEAL,"s","--"),("num_phi",VIOLET,"^","-.")):
        sub=[r for r in conv if r["scan"]==scan]; axes[1].loglog(values(sub,"value"),values(sub,"p95_abs_log10"),color=color,marker=marker,ls=line,lw=1.2,label=scan); axes[2].loglog(values(sub,"value"),values(sub,"runtime_s"),color=color,marker=marker,ls=line,lw=1.2,label=scan)
    panel(axes[1],"b"); finish(axes[1],"grid cells",r"p95 $|\Delta\log F|$","Geometry convergence"); axes[1].legend(fontsize=4.8)
    panel(axes[2],"c"); finish(axes[2],"grid cells","wall time (s)","Measured projection cost")
    save_pub(fig,"fig14_eats_projection")


def fig15_sky_maps() -> None:
    data=np.load(BENCH/"skymap"/"skymap_centroid_motion_data.npz"); summary=rows(BENCH/"skymap"/"skymap_centroid_motion_summary.csv"); decay=np.load(BENCH/"theta_decay"/"theta_j_multiples_bdecay_compare_data.npz")
    fig,axes=plt.subplots(2,4,figsize=(7.15,4.0)); fig.subplots_adjust(left=.06,right=.99,bottom=.17,top=.91,wspace=.38,hspace=.42)
    times=data["skymap_times_s"]; images=[data[f"{dim}_skymap_image_flux"][index] for dim in ("1d","2d") for index in range(len(times))]; positive=np.concatenate([image[image>0] for image in images]); norm=LogNorm(vmin=positive.min(),vmax=max(image.max() for image in images)); mapper=None
    for row_index,dim in enumerate(("1d","2d")):
        x=data[f"{dim}_skymap_x_mas"]; y=data[f"{dim}_skymap_y_mas"]; xx,yy=np.meshgrid(x,y)
        for time_index,time_s in enumerate(times):
            image=data[f"{dim}_skymap_image_flux"][time_index]; total=float(np.sum(image)); xc=float(np.sum(image*xx)/total); yc=float(np.sum(image*yy)/total); sx=float(np.sqrt(np.sum(image*(xx-xc)**2)/total)); sy=float(np.sqrt(np.sum(image*(yy-yc)**2)/total)); ax=axes[row_index,time_index]
            mapper=ax.imshow(image,origin="lower",cmap="magma",norm=norm,extent=[x.min(),x.max(),y.min(),y.max()]); ax.set_xlim(xc-6*sx,xc+6*sx); ax.set_ylim(yc-6*sy,yc+6*sy); ax.plot(xc,yc,marker="+",color="white",ms=5,mew=.8); ax.set_title(rf"{dim.upper()}, $t={time_s:.0e}$ s",fontsize=6.5); ax.set_xlabel(r"$x$ (mas)")
            if time_index==0: ax.set_ylabel(r"$y$ (mas)")
            ax.text(.04,.95,chr(ord("a")+row_index*4+time_index),transform=ax.transAxes,color="white",fontsize=9,fontweight="bold",ha="left",va="top")
    fig.colorbar(mapper,ax=list(axes[:,:3].flat),location="bottom",fraction=.035,pad=.10,label=r"$I_\nu$ (shared absolute scale)")
    for model,color,line in (("1d",BLUE,"-"),("2d",RED,"--")):
        sub=[row for row in summary if row["model"]==model and float(row["theta_v_over_theta_j"])==2 and row["centroid_offset_mas"]]; axes[0,3].loglog(values(sub,"time_s"),values(sub,"centroid_offset_mas"),color=color,ls=line,lw=1.4,label=model.upper())
        speed=[row for row in summary if row["model"]==model and float(row["theta_v_over_theta_j"])==2 and row["beta_app_at_midpoint"]]; axes[1,3].semilogx(values(speed,"time_s"),values(speed,"beta_app_at_midpoint"),color=color,ls=line,lw=1.4,label=model.upper())
    panel(axes[0,3],"d"); finish(axes[0,3],r"$t_{\rm obs}$ (s)","centroid offset (mas)","Centroid motion"); axes[0,3].legend()
    panel(axes[1,3],"h"); finish(axes[1,3],r"$t_{\rm obs}$ (s)",r"$\beta_{\rm app}$","Apparent speed"); axes[1,3].legend(loc="lower left")
    frequency_index=int(np.argmin(np.abs(decay["lc_freqs_hz"]-1e14))); inset=axes[1,3].inset_axes((.48,.52,.50,.43)); inset.loglog(decay["lc_times_s"],decay["flux_1d_lc"][0,frequency_index],color=BLUE,lw=.7); inset.loglog(decay["lc_times_s"],decay["flux_2d_lc"][0,frequency_index],color=RED,ls="--",lw=.7); inset.set_title(r"decay $\alpha_t=-0.4$",fontsize=4.5); inset.tick_params(labelsize=3.8); inset.grid(color=LIGHT,lw=.3)
    save_pub(fig,"fig15_sky_maps")


def fig16_error_budget() -> None:
    table=rows(BENCH/"convergence"/"convergence.csv"); dims=("radius","electron","photon","angle")
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.8)); fig.subplots_adjust(left=.09,right=.98,bottom=.23,top=.84,wspace=.38)
    for index,(dim,color) in enumerate(zip(dims,COLORS)):
        sub=[r for r in table if r["dimension"]==dim and r["component"]=="total"]
        reverse=[r for r in table if r["dimension"]==dim and r["component"]=="reverse"]
        axes[0].loglog(values(sub,"level"),values(sub,"relative_error_p95"),color=color,marker="o",lw=1.2,label=dim)
        axes[0].axhline(float(sub[0]["reference_uncertainty_p95"]),color=color,lw=.55,ls=":",alpha=.65)
        axes[0].loglog(values(reverse,"level"),values(reverse,"relative_error_p95"),color=color,marker="^",lw=.75,ls="--",alpha=.55,label="RS component" if index==0 else None)
        axes[1].loglog(values(sub,"wall_time_median_s"),values(sub,"relative_error_p95"),color=color,marker="o",lw=1.2)
    panel(axes[0],"a"); finish(axes[0],"grid cells","p95 relative flux error","Resolution convergence")
    dimension_legend=axes[0].legend(fontsize=5,ncol=2,loc="upper right"); axes[0].add_artist(dimension_legend)
    axes[0].legend(handles=[Line2D([],[],color=GRAY,marker="o",lw=1.2,label="total"),Line2D([],[],color=GRAY,marker="^",ls="--",lw=.8,label="RS component"),Line2D([],[],color=GRAY,ls=":",lw=.8,label="reference uncertainty")],fontsize=4.7,loc="lower left")
    panel(axes[1],"b"); finish(axes[1],"wall time (s)","p95 relative flux error","Accuracy--cost boundary")
    save_pub(fig,"fig16_error_budget")


def fig17_runtime_reuse() -> None:
    table=rows(BENCH/"runtime"/"runtime_summary.csv"); cases=sorted({(r["medium"],r["dimension"]) for r in table if (r["medium"],r["dimension"]) != ("wind","2d")})
    fig,axes=plt.subplots(2,2,figsize=(7.15,4.2)); fig.subplots_adjust(left=.08,right=.99,bottom=.13,top=.91,wspace=.42,hspace=.52)
    for (case,color) in zip(cases,COLORS):
        sub=sorted([r for r in table if (r["medium"],r["dimension"])==case],key=lambda r:int(r["threads"])); threads=values(sub,"threads")
        label=" ".join(case)
        axes[0,0].loglog(threads,values(sub,"cold_solve_median"),color=color,marker="o",lw=1.2,label=label)
        axes[0,0].loglog(threads,values(sub,"cold_solve_p95"),color=color,lw=.8,ls=":",alpha=.75)
        axes[0,1].loglog(threads,values(sub,"projection_median"),color=color,marker="s",lw=1.2)
        axes[0,1].loglog(threads,values(sub,"projection_p95"),color=color,lw=.8,ls=":",alpha=.75)
        axes[0,1].loglog(threads,values(sub,"warm_query_median"),color=color,marker="^",lw=.9,ls="--",alpha=.75)
        axes[1,1].plot(threads,values(sub,"peak_rss_mib_median"),color=color,marker="o",lw=1.2,label=label)
        axes[1,1].plot(threads,values(sub,"peak_rss_mib_p95"),color=color,lw=.8,ls=":",alpha=.75)

    stages=("cold_solve","projection","warm_query")
    x=np.arange(len(stages))
    for index,stage in enumerate(stages):
        relative=np.asarray([float(row[f"{stage}_iqr"])/float(row[f"{stage}_median"]) for row in table if (row["medium"],row["dimension"]) in cases])
        center=float(np.median(relative)); low=float(np.percentile(relative,25)); high=float(np.percentile(relative,75))
        axes[1,0].errorbar(index,center,yerr=[[center-low],[high-center]],fmt="o",color=COLORS[index],capsize=4,lw=1.2)
    axes[1,0].set_xticks(x,("cold","projection","warm cache"),rotation=15)

    panel(axes[0,0],"a"); finish(axes[0,0],"OpenMP threads","wall time (s)","Cold solve: median and p95"); axes[0,0].legend(fontsize=5)
    panel(axes[0,1],"b"); finish(axes[0,1],"OpenMP threads","wall time (s)","Projection and cached query")
    axes[0,1].legend(handles=[Line2D([],[],color=GRAY,marker="s",lw=1.2,label="projection median"),Line2D([],[],color=GRAY,ls=":",lw=.8,label="projection p95"),Line2D([],[],color=GRAY,marker="^",ls="--",lw=.9,label="warm cached median")],fontsize=4.7)
    panel(axes[1,0],"c"); finish(axes[1,0],"","IQR / median","Run-to-run dispersion")
    panel(axes[1,1],"d"); finish(axes[1,1],"OpenMP threads","peak RSS (MiB)","Memory: median and p95")
    axes[1,1].legend(handles=[Line2D([],[],color=GRAY,marker="o",lw=1.2,label="median"),Line2D([],[],color=GRAY,ls=":",lw=.8,label="p95")],fontsize=4.8)
    save_pub(fig,"fig17_runtime_reuse")


BUILDERS = {
    "fig1": fig1_visual_overview, "fig2": fig2_external_media, "fig3": fig3_radius_state,
    "fig4": fig4_jet_profiles, "fig5": fig5_electron_transport, "fig6": fig6_electron_solvers,
    "fig7": fig7_radiation, "fig8": fig8_primary_reverse_shock, "fig9": fig9_magnetized_reverse_shock,
    "fig10": fig10_wind_termination_shell, "fig11": fig11_density_jumps, "fig12": fig12_smooth_clumpy_medium,
    "fig13": fig13_hadronic_pair_state, "fig14": fig14_eats_projection, "fig15": fig15_sky_maps,
    "fig16": fig16_error_budget, "fig17": fig17_runtime_reuse,
}


def main() -> None:
    parser=argparse.ArgumentParser(description="Render the 17 ASGARD manuscript figures from tracked source data.")
    parser.add_argument("--only",nargs="+",choices=tuple(BUILDERS))
    args=parser.parse_args()
    for name in args.only or BUILDERS: BUILDERS[name]()


if __name__ == "__main__":
    main()
