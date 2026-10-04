"""Distributional diagnostics of the corrected brain-age gap.

The question behind this is whether any single participant is extreme enough to
move a regression coefficient on their own, and whether the two modalities are
equally prone to that. The gap is an out-of-sample prediction error, so nothing
constrains its range: a participant the model has no good analogue for can land
twenty years out without any data problem being present.

z is computed within modality, so |z| > 4 means "four SDs from the mean of that
modality's own gap distribution", not a comparison between modalities. Under a
normal distribution with n = 1,380 the expected count above |z| = 4 is about
0.09, so even one case is more than normality predicts.

Nothing here is a screening rule. No participant is excluded on the basis of
these numbers; they are reported so that a reader can see the tails rather than
infer them from an SD.

OUTPUT  results/table_bag_extremes.csv
        console summary
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results")

RUNS = {
    "Structural":          "BAG_structural_corrected",
    "Functional (PCA)":    "BAG_functional_pca_corrected",
    "Functional (no PCA)": "BAG_functional_nopca_corrected",
}
Z_CUT = 4.0


def main():
    s = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
    rows, z_store = [], {}

    for label, col in RUNS.items():
        g = s[col].to_numpy(float)
        ok = np.isfinite(g)
        g = g[ok]
        m, sd = g.mean(), g.std(ddof=1)
        z = (g - m) / sd
        z_store[label] = pd.Series(z, index=s.loc[ok, "participant_id"].values)

        q = np.percentile(g, [1, 5, 95, 99])
        rows.append({
            "Modality": label,
            "n": int(ok.sum()),
            "M": round(float(m), 2),
            "SD": round(float(sd), 2),
            "Min": round(float(g.min()), 2),
            "P1": round(float(q[0]), 2),
            "P5": round(float(q[1]), 2),
            "P95": round(float(q[2]), 2),
            "P99": round(float(q[3]), 2),
            "Max": round(float(g.max()), 2),
            "Range": round(float(g.max() - g.min()), 2),
            f"n |z|>{Z_CUT:.0f}": int((np.abs(z) > Z_CUT).sum()),
            f"n z<-{Z_CUT:.0f}": int((z < -Z_CUT).sum()),
            f"n z>+{Z_CUT:.0f}": int((z > Z_CUT).sum()),
            "Skew": round(float(pd.Series(g).skew()), 2),
            "Kurtosis (excess)": round(float(pd.Series(g).kurtosis()), 2),
        })

    tbl = pd.DataFrame(rows)
    pd.set_option("display.width", 220)
    print("CORRECTED BRAIN-AGE GAP: RANGE, TAILS AND EXTREME CASES\n")
    print(tbl.to_string(index=False))

    exp = len(z_store["Structural"]) * 2 * (1 - 0.9999683)
    print(f"\nExpected count beyond |z| > {Z_CUT:.0f} under normality: "
          f"{exp:.2f} per modality (n = {len(z_store['Structural'])}).")

    # Are the extreme participants the same people across modalities?
    print("\nOverlap of extreme cases between modalities:")
    labels = list(RUNS)
    sets = {k: set(v.index[np.abs(v) > Z_CUT]) for k, v in z_store.items()}
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            shared = sets[a] & sets[b]
            print(f"  {a} ({len(sets[a])}) vs {b} ({len(sets[b])}): "
                  f"{len(shared)} in common")

    # The same question at a threshold loose enough to have cases everywhere.
    print("\nAt |z| > 3, for context:")
    sets3 = {k: set(v.index[np.abs(v) > 3]) for k, v in z_store.items()}
    for k in labels:
        print(f"  {k}: {len(sets3[k])}")
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            a, b = labels[i], labels[j]
            print(f"  {a} ∩ {b}: {len(sets3[a] & sets3[b])}")

    p = os.path.join(RES, "table_bag_extremes.csv")
    tbl.to_csv(p, index=False)
    print("\nwrote " + os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
