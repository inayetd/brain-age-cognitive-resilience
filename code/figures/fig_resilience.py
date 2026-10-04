"""What cognitive resilience is, and how it is computed.

This is the concept the thesis is named after, so it gets its own figure
and is explained twice: once as a picture and once as a recipe.

  A  the picture. Brain-age gap against resilience. The vertical band
     holds every participant whose gap is between +6 and +10 years - the
     same brain age, in any practical sense - and their cognition spans
     more than four standard deviations. That is the thesis title drawn
     from the data rather than asserted.

  B  the recipe. The model written as a template a non-specialist can
     read, with each term translated underneath, and the leftover named.

WHY THE SCATTER LOOKS UNCORRELATED
  It is, exactly: r = .000. Resilience is defined as the part of
  cognition the gap and the demographics do not explain, so by
  construction it carries none of the gap. That is a property of the
  definition, not a finding, and the figure says so.

WHAT IS NOT CLAIMED
  A person high on this axis is not "protected". They performed better
  than the model expected. The model is cross-sectional and explains
  about 55% of cognitive variance; the rest is measurement, unmeasured
  factors, and chance.

INPUT   results/sem/h3_resilience_scores_<TAG>.csv
        results/bag_sensitivity.csv
OUTPUT  results/figures/fig5a_resilience_concept_{thesis.pdf,presentation.png}
        results/figures/fig5b_resilience_recipe_{thesis.pdf,presentation.png}
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")
SEM = os.path.join(RES, "sem")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))
TAG = os.environ.get("FUNCTIONAL_TAG", "pca_it50")

BLUE = "#15497B"
PURPLE = "#7B2D8E"
GREY = "#5A6570"
LIGHT = "#C7CDD4"
INK = "#1A1A1A"
BAND = (6.0, 10.0)          # the "same brain age" slice, in years
MINUS = "−"

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, lw=1.5, mk=5, alpha=0.26,
                   scale=1.0, title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=15, lw=2.7, mk=13,
                         alpha=0.32, scale=1.5, title_size=19,
                         show_title=True),
}
if os.environ.get("FIG_NO_TITLE"):
    for _s in STYLES.values():
        _s["show_title"] = False


def apply_style(st):
    plt.rcParams.update({
        "figure.dpi": st["dpi"], "font.size": st["base"],
        "axes.labelsize": st["base"], "xtick.labelsize": st["base"] * 0.95,
        "ytick.labelsize": st["base"] * 0.95,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": "#444444",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
        "figure.facecolor": "white",
    })


def load():
    r = pd.read_csv(os.path.join(SEM, f"h3_resilience_scores_{TAG}.csv"))
    b = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    d = (r[r.modality == "structural"]
         .merge(b[["participant_id", "age", "BAG_structural_corrected"]],
                on="participant_id")
         .rename(columns={"BAG_structural_corrected": "gap"}))
    return d.dropna(subset=["gap", "resilience"])


def panel_concept(ax, d, st):
    inside = d[(d.gap > BAND[0]) & (d.gap < BAND[1])]

    ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=0)
    ax.axvspan(BAND[0], BAND[1], color=BLUE, alpha=0.08, zorder=0)

    out = d[(d.gap <= BAND[0]) | (d.gap >= BAND[1])]
    ax.scatter(out.gap, out.resilience, s=st["mk"], alpha=st["alpha"] * 0.65,
               color=GREY, linewidths=0, zorder=2)
    ax.scatter(inside.gap, inside.resilience, s=st["mk"] * 1.25,
               alpha=min(st["alpha"] + 0.25, 0.9), color=BLUE,
               linewidths=0, zorder=3)

    # The two ends of the band, named. These are the thesis in two people.
    hi = inside.loc[inside.resilience.idxmax()]
    lo = inside.loc[inside.resilience.idxmin()]
    for row, txt, va in [(hi, "better cognition\nthan expected", "bottom"),
                         (lo, "worse cognition\nthan expected", "top")]:
        ax.plot([row.gap], [row.resilience], "o", ms=st["mk"] ** 0.5 * 3.1,
                mfc="none", mec=BLUE, mew=st["lw"] * 1.2, zorder=4)
        ax.annotate(txt, xy=(row.gap, row.resilience),
                    xytext=(14, 12 if va == "bottom" else -12),
                    textcoords="offset points", ha="left", va=va,
                    fontsize=st["base"] * 0.86, color=BLUE,
                    linespacing=1.25)

    ax.set_xlabel("Brain-age gap (years)")
    ax.set_ylabel("Cognitive resilience\n(SD of general cognition)")
    ax.set_xlim(d.gap.min() - 2, d.gap.max() + 2)

    ax.annotate(f"{len(inside)} people with a gap of "
                f"+{BAND[0]:.0f} to +{BAND[1]:.0f} years\n"
                f"their cognition spans "
                f"{inside.resilience.max() - inside.resilience.min():.1f} SD",
                xy=(np.mean(BAND), ax.get_ylim()[1]),
                xytext=(0, -6), textcoords="offset points",
                ha="center", va="top", fontsize=st["base"] * 0.9,
                color=BLUE, weight="semibold", linespacing=1.3)

    r = np.corrcoef(d.gap, d.resilience)[0, 1]
    ax.annotate(f"r = {abs(r):.3f} by construction: resilience is what the "
                f"gap does not explain",
                xy=(0.5, -0.185), xycoords="axes fraction", ha="center",
                va="top", fontsize=st["base"] * 0.84, color=GREY)


def panel_recipe(ax, st, hi, lo):
    """The recipe, in three steps, finished with two real participants.

    The worked example is not an illustration: both people are in the data
    and both sit in the +6 to +10 year band of panel A. Showing what the
    number means for an actual pair is what makes the definition land."""
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    B = st["base"]
    W = 96          # not 100: a box flush with the axis has its rounded
                    # corners clipped and renders as two bare lines.

    def box(x, y, w, h, fc, ec, lw=1.0):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.6,rounding_size=1.6",
                                    fc=fc, ec=ec, lw=st["lw"] * lw, zorder=1))

    # --- step 1 -------------------------------------------------------
    ax.text(2, 100, "1.   Predict each person\u2019s cognition from their "
            "brain and their background",
            ha="left", va="top", fontsize=B * 0.95, color=INK,
            weight="semibold")

    box(50, 86, W, 12, "#EFF2F5", "#EFF2F5")
    ax.text(50, 89.0,
            "cognition   ~   brain-age gap  +  age  +  sex  +  education",
            ha="center", va="center", fontsize=B * 1.0, color=INK,
            family="monospace")
    ax.text(50, 82.5,
            "the thing being explained                       what it is "
            "explained by",
            ha="center", va="center", fontsize=B * 0.74, color=GREY)

    ax.add_patch(FancyArrowPatch((50, 79.0), (50, 74.0), arrowstyle="-|>",
                                 mutation_scale=13, color=GREY,
                                 lw=st["lw"], shrinkA=0, shrinkB=0))

    # --- step 2 -------------------------------------------------------
    ax.text(2, 71, "2.   Compare what they actually scored with that "
            "prediction",
            ha="left", va="top", fontsize=B * 0.95, color=INK,
            weight="semibold")

    box(50, 57, W, 14, "white", BLUE, 1.2)
    ax.text(50, 61.5, "what they scored   \u2212   what the model predicted",
            ha="center", va="center", fontsize=B * 0.95, color=INK)
    ax.text(50, 53.5, "=   COGNITIVE RESILIENCE", ha="center", va="center",
            fontsize=B * 1.12, color=BLUE, weight="semibold")

    # --- step 3: two real people -------------------------------------
    ax.text(2, 44, "3.   Two people from this study, with almost the same "
            "brain-age gap",
            ha="left", va="top", fontsize=B * 0.95, color=INK,
            weight="semibold")

    for x, row, word, col, tail in [
            (26, hi, "ABOVE", BLUE, "scored better than predicted"),
            (74, lo, "BELOW", PURPLE, "scored worse than predicted")]:
        box(x, 22, 44, 25, "white", col)
        ax.text(x, 31.0,
                f"age {row.age:.0f}   \u00b7   gap {row.gap:+.1f} years",
                ha="center", va="center", fontsize=B * 0.86, color=GREY)
        ax.text(x, 24.0, f"{row.resilience:+.1f} SD".replace("-", MINUS),
                ha="center", va="center", fontsize=B * 1.6, color=col,
                weight="semibold")
        ax.text(x, 17.0, word + " the prediction", ha="center", va="center",
                fontsize=B * 0.9, color=col, weight="semibold")
        ax.text(x, 12.5, tail, ha="center", va="center",
                fontsize=B * 0.84, color=GREY)

    ax.text(50, 3.5,
            "Their brains look almost equally old. Their cognition does not.\n"
            "That difference is what this thesis measures.",
            ha="center", va="center", fontsize=B * 0.9, color=INK,
            linespacing=1.4)


def main():
    d = load()
    os.makedirs(OUTD, exist_ok=True)
    inside = d[(d.gap > BAND[0]) & (d.gap < BAND[1])]
    hi = inside.loc[inside.resilience.idxmax()]
    lo = inside.loc[inside.resilience.idxmin()]

    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)

        # ---- 5a: the concept, on its own, with room to breathe --------
        fig, ax = plt.subplots(figsize=(8.8 * st["scale"], 6.0 * st["scale"]))
        panel_concept(ax, d, st)
        if st["show_title"]:
            fig.suptitle("Same brain age, different cognition",
                         fontsize=st["title_size"], y=1.0, weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.94])
        else:
            fig.tight_layout()
        p = os.path.join(OUTD, f"fig5a_resilience_concept_{tag}.{st['ext']}")
        fig.savefig(p); plt.close(fig)
        print("wrote " + os.path.relpath(p, ROOT))

        # ---- 5b: the recipe, on its own ------------------------------
        fig, ax = plt.subplots(figsize=(9.6 * st["scale"], 5.8 * st["scale"]))
        panel_recipe(ax, st, hi, lo)
        if st["show_title"]:
            fig.suptitle("How cognitive resilience is computed",
                         fontsize=st["title_size"], y=1.0, weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.94])
        else:
            fig.tight_layout()
        p = os.path.join(OUTD, f"fig5b_resilience_recipe_{tag}.{st['ext']}")
        fig.savefig(p); plt.close(fig)
        print("wrote " + os.path.relpath(p, ROOT))

    print(f"\nband +{BAND[0]:.0f} to +{BAND[1]:.0f} yr: n = {len(inside)}, "
          f"resilience {inside.resilience.min():+.2f} to "
          f"{inside.resilience.max():+.2f}")
    print(f"  high example: age {hi.age:.0f}, gap {hi.gap:+.1f} yr, "
          f"resilience {hi.resilience:+.2f}")
    print(f"  low  example: age {lo.age:.0f}, gap {lo.gap:+.1f} yr, "
          f"resilience {lo.resilience:+.2f}")
    print(f"r(gap, resilience) = {np.corrcoef(d.gap, d.resilience)[0, 1]:+.4f}")


if __name__ == "__main__":
    main()
