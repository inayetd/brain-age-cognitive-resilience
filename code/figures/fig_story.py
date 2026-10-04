"""Figures 1, 2 and 4 of the story sequence.

The sequence teaches the study in order, and each figure answers exactly
one question:

  1  Can brain age be estimated at all, and what did the correction do?
  2  Does the gap matter for cognition?
  3  What is associated with the gap and with resilience?   (fig_h2_h3.py)
  4  What does the whole study amount to?
  5  How much do people differ in the gap?                  (fig_bag_violin.py)

Written for a reader who has never seen the study. That constraint drives
most of the choices here: axes are labelled with what the direction means
rather than with the name of the statistic, the correction is shown as
two numbers and an arrow rather than as a scatter, and nothing is
described in causal language because the design is cross-sectional.

Every number is read from the result files.

OUTPUT  results/figures/fig1_brainage_{thesis.pdf,presentation.png}
        results/figures/fig2_h1_{thesis.pdf,presentation.png}
        results/figures/fig4_summary_{thesis.pdf,presentation.png}
"""
import json
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
COL = {"structural": BLUE, "functional": PURPLE}
NAME = {"structural": "Structural MRI",
        "functional": "Functional connectivity"}
MODS = ["structural", "functional"]
GREY = "#5A6570"
LIGHT = "#C7CDD4"
INK = "#1A1A1A"
MINUS = "−"

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, lw=1.5, ms=6, mk=5,
                   alpha=0.30, scale=1.0, title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=15, lw=2.7, ms=12, mk=13,
                         alpha=0.36, scale=1.5, title_size=19,
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


def mn(x, d=3):
    return f"{x:.{d}f}".replace("-", MINUS)


def save(fig, name, st):
    p = os.path.join(OUTD, f"{name}_{st['tag']}.{st['ext']}")
    fig.savefig(p)
    plt.close(fig)
    print("wrote " + os.path.relpath(p, ROOT))


# --------------------------------------------------------------- the data
def load():
    s = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    met = {}
    for k, p in [("structural", "structural/bag_structural_metrics.json"),
                 ("functional", "functional/bag_metrics_pca.json")]:
        met[k] = json.load(open(os.path.join(RES, p)))
    gap = {"structural": "BAG_structural", "functional": "BAG_functional_pca"}
    h1 = pd.read_csv(os.path.join(SEM, f"h1_estimates_{TAG}.csv"))
    h1 = h1[(h1.spec == "primary") & (h1.type == "standardized") &
            (h1.op == "~")]
    mp = pd.read_csv(os.path.join(SEM, f"multiplicity_{TAG}.csv"))
    return s, met, gap, h1, mp


# ================================================== FIGURE 1  brain age
def figure1(s, met, gap, st):
    # Axes are placed explicitly rather than by gridspec: the two panel
    # headings have to sit in known gaps, and tight_layout moves them.
    fig = plt.figure(figsize=(9.0 * st["scale"], 7.6 * st["scale"]))

    lo = s["age"].min() - 3
    hi = s["age"].max() + 3

    axes = [fig.add_axes([0.09, 0.455, 0.405, 0.42]),
            fig.add_axes([0.545, 0.455, 0.405, 0.42])]

    for j, (ax, k) in enumerate(zip(axes, MODS)):
        pred = s["age"] + s[f"{gap[k]}_raw"]
        ax.plot([lo, hi], [lo, hi], ls="--", color=LIGHT,
                lw=st["lw"] * 0.9, zorder=1)
        ax.scatter(s["age"], pred, s=st["mk"], alpha=st["alpha"],
                   color=COL[k], linewidths=0, zorder=2)
        b = np.polyfit(s["age"], pred, 1)
        ax.plot([lo, hi], np.polyval(b, [lo, hi]), color=COL[k],
                lw=st["lw"] * 1.2, zorder=3)

        m = met[k]
        ax.text(0.045, 0.955,
                f"$R^2$ = {m['outer_r2_mean']:.3f}\n"
                f"MAE = {m['outer_mae_mean']:.2f} years",
                transform=ax.transAxes, ha="left", va="top",
                fontsize=st["base"] * 0.95,
                bbox=dict(boxstyle="round,pad=0.35", fc="white",
                          ec="#BBBBBB", lw=0.8))
        ax.set_title(NAME[k], color=COL[k], fontsize=st["base"] * 1.05,
                     weight="semibold", pad=8)
        ax.set_xlabel("Chronological age (years)")
        if j == 0:
            ax.set_ylabel("Predicted brain age (years)")
        else:
            ax.set_yticklabels([])
            ax.text(0.96, 0.05, "dashed: perfect prediction",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=st["base"] * 0.82, color=GREY)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_aspect("equal", adjustable="box")

    fig.text(0.055, 0.935,
             "A   Structural MRI predicted chronological age more "
             "accurately than functional connectivity",
             ha="left", va="center", fontsize=st["base"] * 1.05,
             weight="semibold")

    # Panel B: the correction as two numbers and an arrow. A scatter would
    # ask the reader to estimate a slope; this asks them only to read.
    fig.text(0.055, 0.315,
             "B   Age-bias correction largely removed the systematic "
             "association between the gap and chronological age",
             ha="left", va="center", fontsize=st["base"] * 1.05,
             weight="semibold")

    ax = fig.add_axes([0.055, 0.03, 0.90, 0.235])
    ax.set_xlim(0, 10); ax.set_ylim(-0.1, 2.5)
    ax.axis("off")

    ax.text(2.75, 2.15, "Before correction", ha="center", va="center",
            fontsize=st["base"] * 0.98, color=GREY, weight="semibold")
    ax.text(6.35, 2.15, "After correction", ha="center", va="center",
            fontsize=st["base"] * 0.98, color=GREY, weight="semibold")

    for i, k in enumerate(MODS):
        y = 1.30 - i * 0.80
        m = met[k]
        ax.text(1.35, y, NAME[k], ha="right", va="center",
                fontsize=st["base"], color=COL[k], weight="semibold")
        ax.text(2.75, y, mn(m["bag_age_corr_raw"], 2), ha="center",
                va="center", fontsize=st["base"] * 1.25, color=COL[k])
        ax.text(6.35, y, f"+{m['bag_age_corr_corrected']:.2f}", ha="center",
                va="center", fontsize=st["base"] * 1.25, color=COL[k])
        ax.add_patch(FancyArrowPatch((3.45, y), (5.60, y),
                                     arrowstyle="-|>", mutation_scale=16,
                                     color=COL[k], lw=st["lw"] * 1.1,
                                     shrinkA=0, shrinkB=0))
    ax.text(9.9, 0.90,
            "correlation between the brain-age gap\nand chronological age",
            ha="right", va="center", fontsize=st["base"] * 0.88, color=GREY,
            linespacing=1.3)

    if st["show_title"]:
        fig.suptitle("Estimating brain age", fontsize=st["title_size"],
                     y=0.995, weight="semibold")
    return fig


# ================================================== FIGURE 2  H1
def figure2(h1, mp, st):
    """H1 shown as the relationship itself rather than as a coefficient.

    A forest plot of beta is exact but abstract: the axis is an association,
    not a quantity, and a reader who has not met standardised coefficients
    before reads the negative sign as a negative gap. Here participants are
    sorted into fifths by their gap and the average cognitive score of each
    fifth is plotted. The axes are then two things anyone can hold: years of
    gap along the bottom, cognitive performance up the side.

    What is plotted is unadjusted, and the H1 test is adjusted for age, sex
    and education; the two agree closely because the corrected gap is almost
    independent of age by construction (r = .02 and .08). The adjusted
    estimate is annotated on each panel so the figure cannot be mistaken for
    the test.

    The scores are factor scores from the measurement model, used here for
    description only. Every inferential model keeps cognition latent.
    """
    gs = pd.read_csv(os.path.join(ROOT, "data", "clean",
                                  "cognitive_factor_scores_V1.csv"))
    bs = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    m = gs.merge(bs, on="participant_id").dropna(subset=["cog_g_score"])
    m["cz"] = (m.cog_g_score - m.cog_g_score.mean()) / m.cog_g_score.std()

    gapcol = {"structural": "BAG_structural_corrected",
              "functional": "BAG_functional_pca_corrected"}
    NQ = 5

    fig, axes = plt.subplots(1, 2, figsize=(6.4 * st["scale"] * 2,
                                            5.4 * st["scale"]), sharey=True)

    for ax, k in zip(axes, MODS):
        d = m.dropna(subset=[gapcol[k]]).copy()
        d["q"] = pd.qcut(d[gapcol[k]], NQ, labels=False)

        xs, ys, es, labs = [], [], [], []
        for i in range(NQ):
            s_ = d[d.q == i]
            c = s_["cz"]
            xs.append(i)
            ys.append(c.mean())
            es.append(1.96 * c.std(ddof=1) / np.sqrt(len(c)))
            labs.append(f"{s_[gapcol[k]].min():+.0f} to "
                        f"{s_[gapcol[k]].max():+.0f}".replace("-", MINUS))

        ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.9, zorder=0)
        ax.plot(xs, ys, "-", color=COL[k], lw=st["lw"] * 1.3, zorder=2)
        ax.errorbar(xs, ys, yerr=es, fmt="o", color=COL[k],
                    ms=st["ms"] * 1.1, lw=st["lw"] * 1.3, capsize=5,
                    markeredgecolor="white", markeredgewidth=st["lw"] * 0.6,
                    zorder=3)

        ax.set_xticks(xs)
        ax.set_xticklabels(labs, fontsize=st["base"] * 0.88)
        ax.set_xlim(-0.55, NQ - 0.45)
        ax.set_xlabel("Brain-age gap (years)\n"
                      "participants split into fifths, youngest-looking "
                      "brains on the left")
        ax.set_title(NAME[k], color=COL[k], fontsize=st["base"] * 1.12,
                     weight="semibold", pad=10)

        r = h1[(h1.modality == k) & (h1.lhs == "g") & (h1.rhs == "BAG")].iloc[0]
        ax.annotate(f"adjusted \u03b2 = {mn(r['est'])}\n"
                    f"[{mn(r['ci.lower'])}, {mn(r['ci.upper'])}]",
                    xy=(0.03, 0.05), xycoords="axes fraction",
                    ha="left", va="bottom", fontsize=st["base"] * 0.9,
                    color=COL[k], linespacing=1.4,
                    bbox=dict(boxstyle="round,pad=0.35", fc="white",
                              ec="#CCCCCC", lw=0.8))

    axes[0].set_ylabel("Average cognitive performance\n"
                       "(SD, sample mean = 0)")
    axes[0].annotate("better than average", xy=(-0.42, 0.255),
                     ha="left", va="center", fontsize=st["base"] * 0.88,
                     color=GREY)
    axes[0].annotate("worse than average", xy=(-0.42, -0.295),
                     ha="left", va="center", fontsize=st["base"] * 0.88,
                     color=GREY)

    fig.text(0.5, -0.085,
             "Each point is the average of about 270 people, with its 95% "
             "interval. The line falls from left to right in both "
             "modalities:\nthe older a brain looks for its age, the lower "
             "that group\u2019s average cognitive performance. The "
             "annotated \u03b2 is the H1 estimate,\nwhich is adjusted for "
             "chronological age, sex and education and models cognition as "
             "a latent variable.",
             ha="center", va="top", fontsize=st["base"] * 0.86, color=GREY,
             linespacing=1.45)

    if st["show_title"]:
        fig.suptitle("The older a brain looks for its age, the lower the "
                     "average cognitive performance",
                     fontsize=st["title_size"], y=1.02, weight="semibold")
    fig.tight_layout()
    return fig


# ============================================ FIGURE 2b  H1, forest version
def figure2b(h1, mp, st):
    """The same result as a coefficient plot.

    Kept alongside the binned version because the two serve different
    readers. This one is the estimate: exact, conventional, and what a
    methods-literate reader wants to see. The binned version shows the
    relationship the estimate summarises, which is what an audience meeting
    the study for the first time needs.

    The direction cues sit BELOW the axis line rather than inside the plot.
    Inside, an arrow reads as a quantity of gap and invites the reader to
    take a negative coefficient for a negative gap; outside, it reads as
    what it is, an annotation of the axis.
    """
    fig, ax = plt.subplots(figsize=(8.4 * st["scale"], 4.6 * st["scale"]))

    rows = [(k, h1[(h1.modality == k) & (h1.lhs == "g") &
                   (h1.rhs == "BAG")].iloc[0]) for k in MODS]

    for j, (k, r) in enumerate(rows):
        y = -j
        ax.plot([r["ci.lower"], r["ci.upper"]], [y, y], color=COL[k],
                lw=st["lw"] * 2.2, solid_capstyle="round", zorder=2)
        ax.plot([r["est"]], [y], "o", color=COL[k], ms=st["ms"] * 1.15,
                markeredgecolor="white", markeredgewidth=st["lw"] * 0.8,
                zorder=3)
        ax.text(r["ci.lower"], y + 0.28,
                f"\u03b2 = {mn(r['est'])}    [{mn(r['ci.lower'])}, "
                f"{mn(r['ci.upper'])}]",
                ha="left", va="bottom", fontsize=st["base"] * 0.98,
                color=COL[k], weight="semibold")

    ax.axvline(0, color="#333333", lw=st["lw"] * 1.2, zorder=1)
    ax.set_yticks([0, -1])
    ax.set_yticklabels([NAME[k] for k, _ in rows])
    for lbl, (k, _) in zip(ax.get_yticklabels(), rows):
        lbl.set_color(COL[k]); lbl.set_weight("semibold")
    ax.tick_params(axis="y", length=0, pad=12)
    ax.spines["left"].set_visible(False)

    lo = min(r["ci.lower"] for _, r in rows)
    ax.set_xlim(lo - 0.035, 0.075)
    ax.set_ylim(-1.55, 0.70)
    ax.set_xlabel("Standardised BAG\u2013cognition association "
                  "(\u03b2, 97.5% CI)", labelpad=10)

    # Direction cues under the axis, drawn in figure space so they cannot be
    # mistaken for data.
    # xy is the arrow HEAD and xytext the tail, so the left cue points left
    # and the right cue points right.
    for tail, head, txt in [
            (0.52, 0.05, "Larger BAG associated\nwith lower cognition"),
            (0.84, 0.99, "Larger BAG associated\nwith higher cognition")]:
        ax.annotate("", xy=(head, -0.30), xytext=(tail, -0.30),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", color=GREY,
                                    lw=st["lw"] * 0.8))
        ax.annotate(txt, xy=((tail + head) / 2, -0.36),
                    xycoords="axes fraction", ha="center", va="top",
                    fontsize=st["base"] * 0.88, color=GREY, linespacing=1.35)

    ax.annotate("Negative \u03b2 values indicate that a 1-SD larger "
                "brain-age gap is associated with lower general cognition.\n"
                "Both associations passed the Bonferroni-corrected threshold "
                "(\u03b1 = .025).",
                xy=(0.5, -0.58), xycoords="axes fraction", ha="center",
                va="top", fontsize=st["base"] * 0.88, color=GREY,
                linespacing=1.4)

    if st["show_title"]:
        fig.suptitle("Larger brain-age gaps were associated with lower "
                     "general cognition",
                     fontsize=st["title_size"], y=1.12, weight="semibold")
        fig.text(0.5, 1.035,
                 "Adjusted for chronological age, sex and education",
                 ha="center", va="top", fontsize=st["base"] * 0.95,
                 color=GREY)
    return fig


# ================================================== FIGURE 4  summary
def figure4(h1, mp, st):
    """The closing slide, written as discussion rather than as results.

    The numbers have already been shown four times by this point in the
    sequence; repeating them here would make this a sixth results figure.
    What a closing slide has to do instead is say what the results mean,
    and each card is therefore a claim plus its implication.

    No causal language and no directed arrows: the design is
    cross-sectional and the final line says so.
    """
    fig, ax = plt.subplots(figsize=(11.4 * st["scale"], 6.5 * st["scale"]))
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
    B = st["base"]

    def card(x, y, w, h, n, heading, body, col):
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                    boxstyle="round,pad=0.7,rounding_size=1.6",
                                    fc="white", ec=col, lw=st["lw"] * 1.1,
                                    zorder=2))
        ax.add_patch(plt.Circle((x - w / 2 + 4.4, y + h / 2 - 4.6), 2.7,
                                fc=col, ec=col, zorder=3))
        ax.text(x - w / 2 + 4.4, y + h / 2 - 4.6, str(n), ha="center",
                va="center", fontsize=B * 0.86, color="white",
                weight="semibold", zorder=4)
        ax.text(x - w / 2 + 9.4, y + h / 2 - 4.6, heading, ha="left",
                va="center", fontsize=B * 1.0, color=col,
                weight="semibold", zorder=4)
        ax.text(x - w / 2 + 4.4, y + h / 2 - 11.0, body, ha="left", va="top",
                fontsize=B * 0.88, color=INK, linespacing=1.52, zorder=4)

    card(26, 75, 46, 37, 1, "Brain age is not one thing",
         "Structural and functional estimates of the\n"
         "same person agree only modestly. Which\n"
         "modality is used changes who is called\n"
         "brain-older, so \u201cbrain age\u201d should not be\n"
         "spoken of without naming the measurement\n"
         "it came from.",
         BLUE)

    card(74, 75, 46, 37, 2, "The gap informs, it does not explain",
         "A larger gap goes with lower cognition over\n"
         "and above chronological age, and with\n"
         "general ability rather than any particular\n"
         "domain. But the association is modest: the\n"
         "gap is a marker of brain status, not an\n"
         "account of cognition.",
         PURPLE)

    card(26, 32, 46, 37, 3, "What the gap leaves out is the point",
         "People whose brains look equally old differ\n"
         "widely in how well they actually perform.\n"
         "That remaining variation is not noise to be\n"
         "minimised; it is the quantity this thesis\n"
         "treats as cognitive resilience, and it is\n"
         "large.",
         BLUE)

    card(74, 32, 46, 37, 4, "Resilience looked material",
         "Of six lifestyle and psychosocial exposures,\n"
         "only household income was associated with\n"
         "both a smaller gap and better-than-expected\n"
         "cognition. Well-being, connectedness and\n"
         "stress \u2014 the mechanisms usually invoked \u2014\n"
         "were not.",
         PURPLE)

    ax.text(50, 5.0,
            "All data are cross-sectional. Every relationship described here "
            "is an association between concurrently measured variables and "
            "licenses no claim about direction.",
            ha="center", va="bottom", fontsize=B * 0.86, color=GREY)

    if st["show_title"]:
        fig.suptitle("Same brain age, different cognition",
                     fontsize=st["title_size"] * 1.1, y=1.015,
                     weight="semibold")
    return fig


def main():
    s, met, gap, h1, mp = load()
    os.makedirs(OUTD, exist_ok=True)
    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        apply_style(st)
        save(figure1(s, met, gap, st), "fig1_brainage", st)
        save(figure2(h1, mp, st), "fig2_h1", st)
        save(figure2b(h1, mp, st), "fig2b_h1_forest", st)
        save(figure4(h1, mp, st), "fig4_summary", st)


if __name__ == "__main__":
    main()
