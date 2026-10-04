"""H2 and H3 as conventional forest plots.

LAYOUT
  Built to the convention used in meta-analysis and epidemiology: a label
  column on the left, the plot in the middle with a line at no effect, and
  the numbers as TEXT on the right. The text column is the part that makes
  a forest plot readable - without it a reader has to estimate each value
  by eye from its position, which is exactly what the earlier version of
  this figure asked them to do.

  Exposures are grouped under bold sub-headings, each with one row per
  modality, so the eye reads down a short indented list rather than across
  twelve equally weighted rows.

DECISIONS ARE CARRIED BY THE MARKER, NOT THE INTERVAL
  Two of the twenty-four intervals exclude zero without surviving
  Benjamini-Hochberg correction. A reader who scans for "does the bar
  cross the line?" would call those findings. So a filled marker means the
  test survived correction and an open marker means it did not, whatever
  the interval does, and the surviving rows are the only ones set in bold.

  Marker area is NOT scaled by precision here. In meta-analysis the square
  is sized by study weight; these are twelve estimates from one sample, so
  there is no weight to show and a varying size would imply one.

TWO FIGURES, NOT ONE
  H2 is in years of gap, H3 in standardised units. Side by side a reader
  compares bar lengths, which across these panels means nothing.

Every number is read from the result files.

INPUT   results/sem/h2_single_estimates_<TAG>.csv
        results/sem/h3_estimates_<TAG>.csv
        results/sem/multiplicity_<TAG>.csv
OUTPUT  results/figures/fig6a_h2_gap_{thesis.pdf,presentation.png}
        results/figures/fig6b_h3_resilience_{thesis.pdf,presentation.png}
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
SEM = os.path.join(RES, "sem")
OUTD = os.environ.get("FIG_OUT", os.path.join(RES, "figures"))
TAG = os.environ.get("FUNCTIONAL_TAG", "pca_it50")

COL = {"structural": "#15497B", "functional": "#7B2D8E"}
MARK = {"structural": "o", "functional": "^"}
SHORT = {"structural": "Structural MRI", "functional": "Functional connectivity"}
MODS = ["structural", "functional"]
GREY = "#5A6570"
RULE = "#B9C0C7"
INK = "#1A1A1A"
Q = 0.05
MINUS = "−"

EXPO = [("income", "Household income"),
        ("sleep", "Sleep quality"),
        ("activity", "Physical activity (high vs low)"),
        ("pwb", "Psychological well-being"),
        ("soc", "Social connectedness"),
        ("stress", "Perceived stress")]

# Layout in a 0-100 coordinate system, so the three columns can be placed
# exactly rather than left to the axes to negotiate.
X_LABEL = 1.0
X_PLOT0, X_PLOT1 = 31.0, 66.0
X_EST = 69.0
X_Q = 97.0

STYLES = {
    "thesis": dict(ext="pdf", dpi=300, base=8.6, lw=1.4, ms=6, scale=1.0,
                   title_size=0, show_title=False),
    "presentation": dict(ext="png", dpi=200, base=14, lw=2.4, ms=10,
                         scale=1.5, title_size=19, show_title=True),
}
if os.environ.get("FIG_NO_TITLE"):
    for _s in STYLES.values():
        _s["show_title"] = False


def bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = np.empty(m)
    q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)


def load():
    mp = pd.read_csv(os.path.join(SEM, f"multiplicity_{TAG}.csv"))
    h2 = pd.read_csv(os.path.join(SEM, f"h2_single_estimates_{TAG}.csv"))
    h3 = pd.read_csv(os.path.join(SEM, f"h3_estimates_{TAG}.csv"))
    h2 = h2[(h2.is_test == True) & (h2.type == "unstandardized_years")]
    h3 = h3[(h3.is_test == True) & (h3.type == "standardized")]

    def attach(d, fam):
        q = mp[mp.family == fam][["modality", "label", "p_adjusted"]]
        d = d.merge(q.rename(columns={"label": "exposure",
                                      "p_adjusted": "q"}),
                    on=["modality", "exposure"], how="left")
        if d["q"].isna().any():
            d["q"] = bh(d.pvalue.values)
            print(f"  NOTE: {fam} not covered by multiplicity file; "
                  "q computed here")
        return d
    return attach(h2, "H2"), attach(h3, "H3")


def fmt(x, d):
    return f"{x:+.{d}f}".replace("-", MINUS)


def qtxt(q):
    return "q < .001" if q < .001 else ("q = " + f"{q:.3f}".replace("0.", "."))


def forest(d, st, dec, xlab, left_txt, right_txt, title, sub, outname):
    B = st["base"]
    # Row plan: a bold sub-heading per exposure, then one row per modality.
    rows = []
    for key, lab in EXPO:
        rows.append(("head", lab, None))
        for m in MODS:
            rows.append(("data", m, d[(d.modality == m) &
                                      (d.exposure == key)].iloc[0]))

    n = len(rows)
    fig, ax = plt.subplots(figsize=(10.6 * st["scale"],
                                    (0.40 * n + 2.5) * st["scale"]))
    ax.set_xlim(0, 100)
    ax.set_ylim(-n - 5.2, 2.3)
    ax.axis("off")

    lo = min(d["ci.lower"].min(), 0)
    hi = max(d["ci.upper"].max(), 0)
    pad = 0.08 * (hi - lo)
    lo, hi = lo - pad, hi + pad
    to_x = lambda v: X_PLOT0 + (v - lo) / (hi - lo) * (X_PLOT1 - X_PLOT0)
    zero = to_x(0)

    # --- column headers -------------------------------------------------
    ax.text(X_LABEL, 2.0, "Exposure", ha="left", va="center",
            fontsize=B * 0.94, color=INK, weight="semibold")
    ax.text(X_EST, 2.0, dec["esthead"], ha="left", va="center",
            fontsize=B * 0.94, color=INK, weight="semibold")
    ax.text(X_Q, 2.0, "q", ha="right", va="center", fontsize=B * 0.94,
            color=INK, weight="semibold", style="italic")
    ax.plot([X_LABEL, X_Q], [1.35, 1.35], color=INK, lw=st["lw"] * 0.8)

    # --- the no-effect line, drawn the full height of the plot ----------
    ax.plot([zero, zero], [-n - 0.1, 0.75], color="#555555",
            lw=st["lw"] * 1.1, zorder=1)

    y = 0.45
    for kind, a, r in rows:
        y -= 1
        if kind == "head":
            ax.text(X_LABEL, y, a, ha="left", va="center",
                    fontsize=B * 0.96, color=INK, weight="semibold")
            continue

        m = a
        sig = r["q"] < Q
        c = COL[m]
        ax.text(X_LABEL + 2.2, y, SHORT[m], ha="left", va="center",
                fontsize=B * 0.90, color=c,
                weight="semibold" if sig else "normal")

        ax.plot([to_x(r["ci.lower"]), to_x(r["ci.upper"])], [y, y], color=c,
                lw=st["lw"] * (1.8 if sig else 1.2),
                alpha=1.0 if sig else 0.70,
                solid_capstyle="butt", zorder=2)
        ax.plot([to_x(r["est"])], [y], MARK[m], color=c if sig else "white",
                ms=st["ms"] * (1.0 if sig else 0.88), markeredgecolor=c,
                markeredgewidth=st["lw"] * (0.6 if sig else 1.15), zorder=3)

        ax.text(X_EST, y,
                f"{fmt(r['est'], dec['d'])}  "
                f"[{fmt(r['ci.lower'], dec['d'])}, "
                f"{fmt(r['ci.upper'], dec['d'])}]",
                ha="left", va="center", fontsize=B * 0.88,
                color=INK if sig else GREY,
                weight="semibold" if sig else "normal")
        ax.text(X_Q, y, qtxt(r["q"]), ha="right", va="center",
                fontsize=B * 0.88, color=c if sig else GREY,
                weight="semibold" if sig else "normal")

    # --- axis under the plot column -------------------------------------
    yax = -n - 0.15
    ticks = np.linspace(lo, hi, 5)
    step = ticks[1] - ticks[0]
    mag = 10 ** np.floor(np.log10(step))
    step = round(step / mag) * mag
    t0 = np.ceil(lo / step) * step
    ticks = np.arange(t0, hi + 1e-9, step)
    ax.plot([to_x(ticks[0]), to_x(ticks[-1])], [yax, yax], color=INK,
            lw=st["lw"] * 0.9)
    for t in ticks:
        ax.plot([to_x(t)] * 2, [yax, yax - 0.14], color=INK,
                lw=st["lw"] * 0.9)
        ax.text(to_x(t), yax - 0.32, f"{t:g}".replace("-", MINUS),
                ha="center", va="top", fontsize=B * 0.86, color=INK)
    ax.text((X_PLOT0 + X_PLOT1) / 2, yax - 0.95, xlab, ha="center", va="top",
            fontsize=B * 0.92, color=INK)

    # --- direction cues, below the axis so they cannot read as data -----
    ax.annotate("", xy=(zero - 13, yax - 1.62), xytext=(zero - 1, yax - 1.62),
                arrowprops=dict(arrowstyle="->", color=GREY,
                                lw=st["lw"] * 0.8))
    ax.text(zero - 7, yax - 1.95, left_txt, ha="center", va="top",
            fontsize=B * 0.84, color=GREY, linespacing=1.3)
    ax.annotate("", xy=(zero + 13, yax - 1.62), xytext=(zero + 1, yax - 1.62),
                arrowprops=dict(arrowstyle="->", color=GREY,
                                lw=st["lw"] * 0.8))
    ax.text(zero + 7, yax - 1.95, right_txt, ha="center", va="top",
            fontsize=B * 0.84, color=GREY, linespacing=1.3)

    legend = ("Filled marker: survives Benjamini-Hochberg correction at "
              "q = .05.   Open marker: does not survive, whatever the "
              "interval does.")
    ax.text(X_LABEL, yax - 3.30, legend, ha="left", va="top",
            fontsize=B * 0.84, color=GREY)
    if sub:
        ax.text(X_LABEL, yax - 3.80, sub, ha="left", va="top",
                fontsize=B * 0.84, color=GREY)

    if st["show_title"]:
        fig.suptitle(title, fontsize=st["title_size"], y=0.995, ha="center",
                     weight="semibold")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
    else:
        fig.tight_layout()

    p = os.path.join(OUTD, f"{outname}_{st['tag']}.{st['ext']}")
    fig.savefig(p)
    plt.close(fig)
    print("wrote " + os.path.relpath(p, ROOT))


def main():
    h2, h3 = load()
    os.makedirs(OUTD, exist_ok=True)

    for tag, st in STYLES.items():
        st = dict(st, tag=tag)
        plt.rcParams.update({"figure.dpi": st["dpi"], "font.size": st["base"],
                             "savefig.bbox": "tight",
                             "savefig.pad_inches": 0.06,
                             "figure.facecolor": "white"})

        forest(h2, st, dict(d=2, esthead="Years of gap per SD  [95% CI]"),
               "Years of brain-age gap per SD of the exposure",
               "smaller,\nyounger-appearing gap",
               "larger,\nolder-appearing gap",
               "Only household income was associated with the brain-age gap",
               "n = 1,380 per modality. Each exposure tested in its own "
               "model, adjusted for age, sex and education.",
               "fig6a_h2_gap")

        forest(h3, st, dict(d=3, esthead="Standardised β  [95% CI]"),
               "Standardised β on latent cognition",
               "worse cognition\nthan expected",
               "better cognition\nthan expected",
               "Only household income was associated with cognitive "
               "resilience",
               "n = 1,342 per modality. Resilience = cognition given the "
               "brain-age gap, age, sex and education. A separate "
               "correction family from H2.",
               "fig6b_h3_resilience")

    print("\nintervals that exclude zero but do not survive correction:")
    for nm, d in [("H2", h2), ("H3", h3)]:
        for _, r in d.iterrows():
            if ((r["ci.lower"] > 0) == (r["ci.upper"] > 0)) and r["q"] >= Q:
                print(f"  {nm} {r['modality']:11s} {r['exposure']:9s} "
                      f"q = {r['q']:.3f}")


if __name__ == "__main__":
    main()
