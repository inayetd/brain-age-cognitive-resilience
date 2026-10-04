"""H1: the brain-age gap and latent cognitive performance.

TWO FIGURES, BECAUSE THEY ARE FOR DIFFERENT PLACES

fig10_h1  The result. Two coefficients, nothing else.
    Age, sex and education are in the model but are NOT in this plot. They
    are adjusters, not results, and putting them on the same axis invites
    the reader to compare the hypothesis against its own controls, which
    is not the comparison the figure is for. One message: a larger gap
    goes with lower general cognition in both modalities, more strongly
    for the structural gap.

fig11_h1_supporting  The backup slide.
    Two things worth having ready when asked, and only then:
    the gap expressed in years of chronological ageing, and the fact that
    it tracks general ability rather than executive- or memory-specific
    variance. Neither is the H1 result and neither belongs on the results
    slide.

ON THE INTERVALS
  H1 was tested at alpha = .025 with 97.5% intervals, Bonferroni across
  the two modalities. The intervals drawn are those, and every axis label
  says so.

INPUT   results/sem/h1_estimates_<TAG>.csv
        results/bag_sensitivity.csv      (SDs, for the year conversion)
OUTPUT  results/figures/fig10_h1_{thesis.pdf,presentation.png}
        results/figures/fig11_h1_supporting_{thesis.pdf,presentation.png}
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))
TAG = os.environ.get("FUNCTIONAL_TAG", "pca_it50")

COL = {"structural": "#15497B", "functional": "#7B2D8E"}
LABEL = {"structural": "Structural BAG", "functional": "Functional (PCA) BAG"}
MODS = ["structural", "functional"]
GREY = "#5A6570"
CONF_PCT = 97.5
ALPHA = 0.025

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, lw=1.5, ms=6, scale=1.0,
                   title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=16, lw=2.8, ms=13,
                         scale=1.55, title_size=20, show_title=True),
}


# On a slide the heading is set in the deck's own typography, so the figure
# must not carry one as well. FIG_NO_TITLE=1 suppresses it for that build.
if os.environ.get("FIG_NO_TITLE"):
    for _st in STYLES.values():
        _st["show_title"] = False


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"],
        "font.size": st["base"],
        "axes.labelsize": st["base"],
        "xtick.labelsize": st["base"] * 0.95,
        "ytick.labelsize": st["base"] * 1.05,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "axes.edgecolor": "#444444",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.04,
        "figure.facecolor": "white",
    })


def load():
    f = os.path.join(RES, "sem", f"h1_estimates_{TAG}.csv")
    if not os.path.exists(f):
        raise SystemExit(f"not found: {f}\n  Run 04_h1_base.R first.")
    e = pd.read_csv(f)
    e = e[(e.spec == "primary") & (e.type == "standardized") & (e.op == "~")]
    s = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    sd = {"structural": s["BAG_structural_corrected"].std(ddof=1),
          "functional": s["BAG_functional_pca_corrected"].std(ddof=1)}
    return e, sd, s["age"].std(ddof=1)


def get(e, modality, lhs, rhs):
    r = e[(e.modality == modality) & (e.lhs == lhs) & (e.rhs == rhs)]
    if not len(r):
        raise SystemExit(f"missing coefficient: {modality} {lhs} ~ {rhs}")
    r = r.iloc[0]
    return r["est"], r["ci.lower"], r["ci.upper"], r["pvalue"]


def fmt_p(p):
    return "p < .001" if p < .001 else f"p = {p:.3f}".replace("0.", ".")


def mn(x, d=3):
    """A real minus sign, so the annotation matches the axis ticks."""
    return f"{x:.{d}f}".replace("-", "\u2212")


# ---------------------------------------------------------------------
# fig10: the result
# ---------------------------------------------------------------------
def fig_result(e, st):
    fig, ax = plt.subplots(figsize=(7.4 * st["scale"], 3.5 * st["scale"]))

    for j, m in enumerate(MODS):
        b, lo, hi, p = get(e, m, "g", "BAG")
        y = -j
        ax.plot([lo, hi], [y, y], color=COL[m], lw=st["lw"] * 1.6,
                solid_capstyle="round", zorder=2)
        ax.plot([b], [y], "o", color=COL[m], ms=st["ms"],
                markeredgecolor="white", markeredgewidth=st["lw"] * 0.6,
                zorder=3)
        # Fixed x rather than hi + offset, so the two annotations align in a
        # column instead of stepping with the interval width.
        ax.text(0.022, y,
                f"β = {mn(b)}   [{mn(lo)}, {mn(hi)}]   {fmt_p(p)}",
                va="center", ha="left", fontsize=st["base"] * 0.95,
                color=COL[m])

    ax.axvline(0, color="#333333", lw=st["lw"], zorder=1)
    ax.set_yticks([-j for j in range(len(MODS))])
    ax.set_yticklabels([LABEL[m] for m in MODS])
    for lbl, m in zip(ax.get_yticklabels(), MODS):
        lbl.set_color(COL[m])
        lbl.set_weight("semibold")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel(f"Standardised β on general cognition  "
                  f"({CONF_PCT}% CI)")
    ax.set_xlim(-0.285, 0.40)
    ax.set_ylim(-1.75, 0.75)

    ax.annotate(f"Bonferroni threshold: α = {ALPHA:.3f}".replace("0.", "."),
                xy=(0.012, 0.04), xycoords="axes fraction",
                ha="left", va="bottom",
                fontsize=st["base"] * 0.88, color=GREY)

    if st["show_title"]:
        fig.suptitle("Higher brain-age gap, lower general cognition",
                     fontsize=st["title_size"], y=1.0, ha="center",
                     weight="semibold")
        fig.tight_layout(rect=[0, 0, 1, 0.90])
    else:
        fig.tight_layout()
    return fig


# ---------------------------------------------------------------------
# fig11: the backup slide
# ---------------------------------------------------------------------
def fig_supporting(e, sd, age_sd, st):
    fig, axes = plt.subplots(1, 2, figsize=(5.2 * st["scale"] * 2,
                                            4.0 * st["scale"]))

    # (a) one year of gap, in years of chronological ageing. Both are in
    # years, so the ratio of the per-year effects has no units left in it.
    # The interval comes from the coefficient; the two SDs are treated as
    # known, which understates the uncertainty slightly.
    ax = axes[0]
    vals = {}
    for m in MODS:
        b, lo, hi, _ = get(e, m, "g", "BAG")
        a, _, _, _ = get(e, m, "g", "age_z")
        per_age = a / age_sd
        vals[m] = tuple(x / sd[m] / per_age for x in (b, lo, hi))

    for j, m in enumerate(MODS):
        r, lo, hi = vals[m]
        y = -j
        ax.barh([y], [r], height=0.40, color=COL[m], alpha=0.28,
                edgecolor=COL[m], lw=st["lw"], zorder=2)
        ax.plot([min(lo, hi), max(lo, hi)], [y, y], color=COL[m],
                lw=st["lw"] * 1.4, solid_capstyle="round", zorder=3)
        ax.text(max(r, max(lo, hi)) + 0.04, y,
                f"{r:.2f} yr  [{min(lo, hi):.2f}, {max(lo, hi):.2f}]",
                va="center", ha="left", fontsize=st["base"] * 0.95,
                color=COL[m])
    ax.set_yticks([0, -1])
    ax.set_yticklabels([LABEL[m] for m in MODS])
    for lbl, m in zip(ax.get_yticklabels(), MODS):
        lbl.set_color(COL[m]); lbl.set_weight("semibold")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel("Years of chronological ageing equivalent\n"
                  "to one year of brain-age gap")
    ax.set_xlim(0, max(max(v) for v in vals.values()) * 1.85)
    ax.set_ylim(-1.6, 0.6)
    ax.set_title("a   What one year of gap is worth", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    # (b) the gap against the general factor and the two specific factors.
    ax = axes[1]
    targets = [("g", "General\ncognition"), ("Exec", "Executive\nspecific"),
               ("Mem", "Memory\nspecific")]
    off = 0.17
    for i, (lhs, _) in enumerate(targets):
        for j, m in enumerate(MODS):
            b, lo, hi, _ = get(e, m, lhs, "BAG")
            y = -i + (off if j == 0 else -off)
            ax.plot([lo, hi], [y, y], color=COL[m], lw=st["lw"] * 1.3,
                    solid_capstyle="round", zorder=2)
            # Hollow where the interval covers zero, so the eye does not
            # have to check the axis.
            filled = (lo > 0) == (hi > 0)
            ax.plot([b], [y], "o", color=COL[m] if filled else "white",
                    ms=st["ms"] * 0.85, markeredgecolor=COL[m],
                    markeredgewidth=st["lw"], zorder=3)
    ax.axvline(0, color="#333333", lw=st["lw"] * 0.8, zorder=1)
    ax.set_yticks([-i for i in range(len(targets))])
    ax.set_yticklabels([t[1] for t in targets])
    ax.tick_params(axis="y", length=0, pad=8)
    ax.set_xlabel(f"Standardised β of the gap  ({CONF_PCT}% CI)")
    ax.set_ylim(-len(targets) + 0.5, 0.5)
    ax.annotate("filled: interval excludes zero",
                xy=(0.98, 0.04), xycoords="axes fraction",
                ha="right", va="bottom",
                fontsize=st["base"] * 0.85, color=GREY)
    ax.set_title("b   General ability, not a specific domain", loc="left",
                 fontsize=st["base"] * 1.02, weight="semibold")

    for j, m in enumerate(MODS):
        fig.text(0.015 + j * 0.16, -0.03, LABEL[m], color=COL[m],
                 fontsize=st["base"], weight="semibold")

    if st["show_title"]:
        fig.suptitle("Supporting detail for H1",
                     fontsize=st["title_size"], y=1.0, ha="center",
                     weight="semibold")
        fig.tight_layout(rect=[0, 0, 1, 0.93])
    else:
        fig.tight_layout()
    return fig


def main():
    e, sd, age_sd = load()
    os.makedirs(OUTD, exist_ok=True)

    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        for name, fig in [("fig10_h1", fig_result(e, st)),
                          ("fig11_h1_supporting",
                           fig_supporting(e, sd, age_sd, st))]:
            p = os.path.join(OUTD, f"{name}_{tag}.{st['ext']}")
            fig.savefig(p)
            plt.close(fig)
            print("wrote " + os.path.relpath(p, ROOT))

    print("\nnumbers for the slide:")
    for m in MODS:
        b, lo, hi, p = get(e, m, "g", "BAG")
        print(f"  {LABEL[m]:24s} beta = {b:.3f}  "
              f"[{lo:.3f}, {hi:.3f}]  {fmt_p(p)}")
    print("\nyears of chronological ageing per year of gap:")
    for m in MODS:
        b, lo, hi, _ = get(e, m, "g", "BAG")
        a, _, _, _ = get(e, m, "g", "age_z")
        pa = a / age_sd
        r, l, h = (x / sd[m] / pa for x in (b, lo, hi))
        print(f"  {m:11s} {r:.2f}  [{min(l, h):.2f}, {max(l, h):.2f}]")


if __name__ == "__main__":
    main()
