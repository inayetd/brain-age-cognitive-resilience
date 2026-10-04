"""The income association as distributions, grouped violins.

WHY THIS FIGURE NEEDS A WARNING ON IT
  The association is real and survives correction, but it is small next to
  the spread of the gap: about one year of difference between the lowest
  and highest income third, against a between-person SD of 7.4 years. Drawn
  as violins alone, the three groups would look identical and the figure
  would deny a result the data support.

  So the violins show the distributions, which do overlap almost entirely,
  and the group mean with its 95% interval is drawn on top, which is where
  the shift lives. Both are true and a reader needs both: the effect is
  consistent and it is small.

WHAT THIS IS NOT
  The confirmatory tests treated income as continuous. Splitting it into
  thirds is for display only and throws information away; the coefficient
  from the actual model is annotated so the figure cannot be read as the
  test.

INPUT   results/bag_sensitivity.csv
        data/clean/analysis_table_pca.csv
        results/sem/h3_resilience_scores_<TAG>.csv
        results/sem/multiplicity_<TAG>.csv
OUTPUT  results/figures/fig7_income_{thesis.pdf,presentation.png}
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")
SEM = os.path.join(RES, "sem")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))
TAG = os.environ.get("FUNCTIONAL_TAG", "pca_it50")

COL = {"structural": "#15497B", "functional": "#7B2D8E"}
NAME = {"structural": "Structural MRI", "functional": "Functional connectivity"}
MODS = ["structural", "functional"]
GREY = "#5A6570"
LIGHT = "#C7CDD4"
INK = "#1A1A1A"
MINUS = "−"
GROUPS = ["Lowest third", "Middle third", "Highest third"]

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=9.0, lw=1.5, ms=6, scale=1.0,
                   title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=15, lw=2.6, ms=11,
                         scale=1.5, title_size=19, show_title=True),
}
if os.environ.get("FIG_NO_TITLE"):
    for _s in STYLES.values():
        _s["show_title"] = False


def load():
    d = pd.read_csv(os.path.join(ROOT, "data", "clean",
                                 "analysis_table_pca.csv"))
    b = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    r = pd.read_csv(os.path.join(SEM, f"h3_resilience_scores_{TAG}.csv"))
    mp = pd.read_csv(os.path.join(SEM, f"multiplicity_{TAG}.csv"))

    m = d[["participant_id", "income_log"]].merge(b, on="participant_id")
    m = m.dropna(subset=["income_log"])
    m["g"] = pd.qcut(m.income_log, 3, labels=GROUPS)

    res = {}
    for k in MODS:
        rr = r[r.modality == k][["participant_id", "resilience"]]
        res[k] = m[["participant_id", "g"]].merge(rr, on="participant_id")
    return m, res, mp


def panel(ax, get_series, st, ylab, coef_txt):
    W = 0.34
    for gi, g in enumerate(GROUPS):
        for mi, k in enumerate(MODS):
            v = get_series(g, k)
            x = gi + (-W / 2 if mi == 0 else W / 2) * 1.05

            parts = ax.violinplot([v], positions=[x], widths=W * 0.95,
                                  showmeans=False, showmedians=False,
                                  showextrema=False)
            for bd in parts["bodies"]:
                bd.set_facecolor(COL[k]); bd.set_alpha(0.32)
                bd.set_edgecolor(COL[k]); bd.set_linewidth(st["lw"] * 0.7)

            q1, med, q3 = np.percentile(v, [25, 50, 75])
            ax.add_patch(plt.Rectangle((x - 0.022, q1), 0.044, q3 - q1,
                                       facecolor="white", edgecolor=COL[k],
                                       lw=st["lw"] * 0.8, zorder=3))
            ax.plot([x - 0.042, x + 0.042], [med, med], color=COL[k],
                    lw=st["lw"] * 1.1, zorder=4, solid_capstyle="butt")

            # The mean with its interval, on top: this is where the shift is.
            mu = v.mean()
            ci = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
            ax.errorbar([x], [mu], yerr=[ci], fmt="D", color=COL[k],
                        ms=st["ms"] * 0.62, lw=st["lw"] * 1.5, capsize=4,
                        markeredgecolor="white",
                        markeredgewidth=st["lw"] * 0.5, zorder=6)
            # Against a spread this wide the shift in the mean cannot be
            # seen, so it is also written down.
            ax.annotate(f"{mu:+.2f}".replace("-", MINUS),
                        xy=(x, mu), xytext=(0, 13), textcoords="offset points",
                        ha="center", va="bottom", fontsize=st["base"] * 0.80,
                        color=COL[k], weight="semibold", zorder=7)

    ax.axhline(0, color=LIGHT, lw=st["lw"] * 0.8, zorder=0)
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GROUPS)
    ax.set_xlim(-0.55, len(GROUPS) - 0.45)
    ax.set_ylabel(ylab)
    ax.annotate(coef_txt, xy=(0.5, -0.155), xycoords="axes fraction",
                ha="center", va="top", fontsize=st["base"] * 0.86,
                color=GREY, linespacing=1.45)


def main():
    m, res, mp = load()
    os.makedirs(OUTD, exist_ok=True)

    def q2(k):
        r = mp[(mp.family == "H2") & (mp.modality == k) &
               (mp.label == "income")].iloc[0]
        return r["estimate"], r["p_adjusted"]

    def q3(k):
        r = mp[(mp.family == "H3") & (mp.modality == k) &
               (mp.label == "income")].iloc[0]
        return r["estimate"], r["p_adjusted"]

    gapcol = {"structural": "BAG_structural_corrected",
              "functional": "BAG_functional_pca_corrected"}

    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        plt.rcParams.update({
            "figure.dpi": st["dpi"], "font.size": st["base"],
            "axes.labelsize": st["base"], "xtick.labelsize": st["base"],
            "ytick.labelsize": st["base"] * 0.95,
            "axes.spines.top": False, "axes.spines.right": False,
            "axes.edgecolor": "#444444", "savefig.bbox": "tight",
            "savefig.pad_inches": 0.05, "figure.facecolor": "white"})

        fig, axes = plt.subplots(1, 2, figsize=(5.3 * st["scale"] * 2,
                                                5.9 * st["scale"]))

        e_s, _ = q2("structural"); e_f, _ = q2("functional")
        panel(axes[0],
              lambda g, k: m.loc[m.g == g, gapcol[k]].dropna().values, st,
              "Corrected brain-age gap (years)",
              f"model estimate per SD of log income:\n"
              f"{e_s:+.2f} yr (structural), {e_f:+.2f} yr (functional)")
        axes[0].set_title("A   Brain-age gap", loc="left",
                          fontsize=st["base"] * 1.1, weight="semibold")

        b_s, _ = q3("structural"); b_f, _ = q3("functional")
        panel(axes[1],
              lambda g, k: res[k].loc[res[k].g == g, "resilience"]
                                 .dropna().values, st,
              "Cognitive resilience (SD)",
              f"model estimate: β = {b_s:+.3f} (structural) and "
              f"{b_f:+.3f} (functional) per SD of log income")
        axes[1].set_title("B   Cognitive resilience", loc="left",
                          fontsize=st["base"] * 1.1, weight="semibold")

        for ax in axes:
            ax.set_xlabel("Household income", labelpad=8)

        handles = [Line2D([], [], color=COL[k], marker="s", ls="none",
                          ms=st["ms"] * 0.8, alpha=0.6, label=NAME[k])
                   for k in MODS]
        handles.append(Line2D([], [], color=INK, marker="D", ls="none",
                              ms=st["ms"] * 0.62,
                              label="group mean with 95% interval"))
        fig.legend(handles=handles, loc="lower center", ncol=3,
                   frameon=False, bbox_to_anchor=(0.5, -0.075),
                   fontsize=st["base"] * 0.92)

        # Wrapped by hand: a single long line forces bbox_inches="tight"
        # to widen the whole figure to contain it.
        fig.text(0.5, -0.125,
                 "The distributions overlap almost entirely. The shift "
                 "between the lowest and highest income third is about one "
                 "year of gap,\n"
                 "against a between-person SD of 7.4 years. The association "
                 "is consistent across modalities and survives correction, "
                 "and it is small.\n"
                 "Income was treated as continuous in the models; the "
                 "thirds are for display only.",
                 ha="center", va="top", fontsize=st["base"] * 0.86,
                 color=GREY, linespacing=1.5)

        if st["show_title"]:
            fig.suptitle("Higher household income, a slightly younger brain "
                         "and slightly better cognition than expected",
                         fontsize=st["title_size"], y=1.0, ha="center",
                         weight="semibold")
            fig.tight_layout(rect=[0, 0, 1, 0.94])
        else:
            fig.tight_layout()

        p = os.path.join(OUTD, f"fig7_income_{tag}.{st['ext']}")
        fig.savefig(p)
        plt.close(fig)
        print("wrote " + os.path.relpath(p, ROOT))

    print("\ngroup means:")
    for k in MODS:
        print(f"  {k}")
        for g in GROUPS:
            v = m.loc[m.g == g, gapcol[k]].dropna()
            w = res[k].loc[res[k].g == g, "resilience"].dropna()
            print(f"    {g:14s} gap {v.mean():+6.2f} yr   "
                  f"resilience {w.mean():+5.2f} SD   n = {len(v)}")


if __name__ == "__main__":
    main()
