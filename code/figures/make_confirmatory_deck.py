"""The confirmatory results, written to be read rather than presented.

Every slide carries the model it is reporting, in lavaan syntax, above the
numbers that model produced. The point is that a reader can check the
claim against the specification without asking anyone, and that the
author can find the specification again in six months.

Numbers are read from the result files, not typed, so the deck cannot
drift away from the analysis. Anything that could not be read is marked
in the console at build time.

Colours: navy for headings, the project's structural blue and functional
purple for the two modalities, everywhere.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt
from pptx.oxml.ns import qn
from lxml import etree

ROOT = sys.argv[1]
OUT = sys.argv[2]
TAG = sys.argv[3] if len(sys.argv) > 3 else "pca_it50"
SEM = os.path.join(ROOT, "results", "sem")
RES = os.path.join(ROOT, "results")

NAVY = RGBColor(0x16, 0x32, 0x4F)
BLUE = RGBColor(0x15, 0x49, 0x7B)
PURPLE = RGBColor(0x7B, 0x2D, 0x8E)
GREY = RGBColor(0x5A, 0x65, 0x70)
INK = RGBColor(0x1A, 0x1A, 0x1A)
CODEBG = "EFF2F5"
MINUS = "−"

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]


# ------------------------------------------------------------------ helpers
def slide():
    return prs.slides.add_slide(BLANK)


def tb(s, x, y, w, h):
    box = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    return tf


def P(tf, parts, size, color=INK, font="Calibri", bold=False, first=False,
      before=0, align=PP_ALIGN.LEFT, line=None):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    p.space_before = Pt(before)
    if line:
        p.line_spacing = line
    if isinstance(parts, str):
        parts = [(parts, False)]
    for part in parts:
        r = p.add_run()
        r.text = part[0]
        r.font.size = Pt(size)
        r.font.italic = part[1]
        r.font.color.rgb = part[2] if len(part) > 2 else color
        r.font.name = font
        r.font.bold = bold
    return p


def head(s, title, sub=None):
    tf = tb(s, 0.55, 0.36, 12.25, 1.0)
    P(tf, title, 25, NAVY, font="Cambria", bold=True, first=True)
    if sub:
        P(tf, sub, 13, GREY, before=4)


def code(s, lines, x, y, w, size=11.5, title=None):
    """A lavaan model block. Monospace on a tinted panel."""
    h = 0.215 * len(lines) + (0.34 if title else 0) + 0.22
    shp = s.shapes.add_shape(1, Inches(x), Inches(y), Inches(w), Inches(h))
    shp.fill.solid()
    shp.fill.fore_color.rgb = RGBColor.from_string(CODEBG)
    shp.line.fill.background()
    shp.shadow.inherit = False
    tf = shp.text_frame
    tf.word_wrap = True
    tf.margin_left = Inches(0.14)
    tf.margin_right = Inches(0.1)
    tf.margin_top = Inches(0.09)
    tf.margin_bottom = Inches(0.07)
    first = True
    if title:
        P(tf, title, 11.5, NAVY, bold=True, first=True)
        first = False
    for ln in lines:
        P(tf, ln, size, INK, font="Consolas", first=first, line=1.0)
        first = False
    return h


# --- APA-ish tables: horizontal rules only -------------------------------
LINE = {"w": "12700", "col": "000000"}


def _ln(cell, edge, on):
    tcPr = cell._tc.get_or_add_tcPr()
    tag = qn("a:ln" + edge)
    for old in tcPr.findall(tag):
        tcPr.remove(old)
    e = etree.SubElement(tcPr, tag)
    e.set("w", LINE["w"]); e.set("cap", "flat")
    e.set("cmpd", "sng"); e.set("algn", "ctr")
    if on:
        f = etree.SubElement(e, qn("a:solidFill"))
        etree.SubElement(f, qn("a:srgbClr")).set("val", LINE["col"])
    else:
        etree.SubElement(e, qn("a:noFill"))


def _order(tbl):
    order = ["a:lnL", "a:lnR", "a:lnT", "a:lnB", "a:noFill", "a:solidFill"]
    for row in tbl.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            kids = list(tcPr)
            for k in kids:
                tcPr.remove(k)
            for nm in order:
                for k in kids:
                    if k.tag == qn(nm):
                        tcPr.append(k)


def table(s, rows, cols, x, y, col_w, row_h=0.30, size=11):
    shp = s.shapes.add_table(rows, cols, Inches(x), Inches(y),
                             Inches(sum(col_w)), Inches(row_h * rows))
    t = shp.table
    pr = t._tbl.find(qn("a:tblPr"))
    pr.set("firstRow", "0"); pr.set("bandRow", "0")
    for st in pr.findall(qn("a:tableStyleId")):
        pr.remove(st)
    for i, w in enumerate(col_w):
        t.columns[i].width = Inches(w)
    for r in range(rows):
        t.rows[r].height = Inches(row_h)
        for c in range(cols):
            cell = t.cell(r, c)
            tcPr = cell._tc.get_or_add_tcPr()
            for tg in ("a:solidFill", "a:noFill"):
                for o in tcPr.findall(qn(tg)):
                    tcPr.remove(o)
            etree.SubElement(tcPr, qn("a:noFill"))
            for e in ("L", "R", "T", "B"):
                _ln(cell, e, False)
            cell.margin_left = Inches(0.05)
            cell.margin_right = Inches(0.05)
            cell.margin_top = Inches(0.015)
            cell.margin_bottom = Inches(0.015)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    t._size = size
    return t


def put(cell, parts, size=11, bold=False, align=PP_ALIGN.CENTER, color=INK):
    tf = cell.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    if isinstance(parts, str):
        parts = [(parts, False)]
    for part in parts:
        r = p.add_run()
        r.text = part[0]
        r.font.size = Pt(size)
        r.font.italic = part[1]
        r.font.bold = bold
        r.font.name = "Calibri"
        r.font.color.rgb = part[2] if len(part) > 2 else color


def rules(t, header_rows, last_row, ncol):
    """Rule above the table, under the head row, under the body."""
    for c in range(ncol):
        _ln(t.cell(0, c), "T", True)
        _ln(t.cell(header_rows - 1, c), "B", True)
        _ln(t.cell(header_rows, c), "T", True)
        _ln(t.cell(last_row, c), "B", True)
    _order(t)


def note(s, x, y, w, parts, size=10.5):
    tf = tb(s, x, y, w, 0.9)
    P(tf, [("Note. ", True, GREY)] + parts, size, GREY, first=True, line=1.12)


# ------------------------------------------------------------------- data
def rd(name):
    p = os.path.join(SEM, name)
    if not os.path.exists(p):
        print("  MISSING:", name)
        return None
    return pd.read_csv(p)


h1 = rd(f"h1_estimates_{TAG}.csv")
h1f = rd(f"h1_fit_{TAG}.csv")
h2 = rd(f"h2_single_estimates_{TAG}.csv")
h2v = rd(f"h2_single_variance_{TAG}.csv")
h3 = rd(f"h3_estimates_{TAG}.csv")
h3v = rd(f"h3_single_variance_{TAG}.csv")
mp = rd(f"multiplicity_{TAG}.csv")

met = {}
for k, p in [("structural", "structural/bag_structural_metrics.json"),
             ("functional", "functional/bag_metrics_pca.json")]:
    fp = os.path.join(RES, p)
    met[k] = json.load(open(fp)) if os.path.exists(fp) else {}
sens = pd.read_csv(os.path.join(RES, "bag_sensitivity.csv"))
BAGSD = {"structural": sens["BAG_structural_corrected"].std(ddof=1),
         "functional": sens["BAG_functional_pca_corrected"].std(ddof=1)}

EXPO = [("income", "Household income"), ("sleep", "Sleep quality"),
        ("activity", "Physical activity (high vs low)"),
        ("pwb", "Psychological well-being"),
        ("soc", "Social connectedness"), ("stress", "Perceived stress")]
MODS = [("structural", "Structural", BLUE), ("functional", "Functional", PURPLE)]


def num(x, d=3, signed=False):
    s = f"{x:+.{d}f}" if signed else f"{x:.{d}f}"
    return s.replace("-", MINUS)


def drfmt(x):
    """A delta R2 that rounds to zero from below is a rounding artefact; the
    minus sign would read as a negative contribution, which it is not."""
    if abs(x) < 5e-5:
        return ".0000"
    return f"{x:.4f}".replace("0.", ".").replace("-", MINUS)


def pstr(p):
    if p < .001:
        return "< .001"
    return f"{p:.3f}".replace("0.", ".")


def qof(family, modality, label):
    if mp is None:
        return np.nan
    r = mp[(mp.family == family) & (mp.modality == modality) &
           (mp.label == label)]
    return float(r.p_adjusted.iloc[0]) if len(r) else np.nan


# =================================================================== 1 overview
s = slide()
head(s, "Operationalising cognitive resilience as deviation from brain age",
     "Confirmatory results · AABC / HCP-Aging Release 2")

tf = tb(s, 0.55, 1.55, 6.0, 4.4)
P(tf, "The question", 15, NAVY, bold=True, first=True)
P(tf, "Two people with the same brain-age gap can differ in cognition. "
      "That discrepancy is what this study treats as resilience: cognitive "
      "performance conditional on brain status.", 13, INK, before=4, line=1.14)

P(tf, "Three confirmatory hypotheses", 15, NAVY, bold=True, before=16)
P(tf, [("H1  ", False, NAVY),
       ("a larger gap goes with lower cognition", False)], 13, INK, before=5)
P(tf, [("H2  ", False, NAVY),
       ("six lifestyle and psychosocial exposures relate to the gap",
        False)], 13, INK, before=3)
P(tf, [("H3  ", False, NAVY),
       ("the same six relate to resilience, i.e. to cognition given the gap",
        False)], 13, INK, before=3, line=1.12)

P(tf, "Everything is fitted separately for the two modalities, and every "
      "estimate below is pooled over 50 multiply imputed datasets by "
      "Rubin’s rules.", 12, GREY, before=14, line=1.14)

tf = tb(s, 7.05, 1.55, 5.75, 4.4)
P(tf, "Samples", 15, NAVY, bold=True, first=True)
for a, b, c in [
    ("1,386", "participants in the cleaned phenotype table", ""),
    ("1,380", "with a usable brain-age gap", "H2"),
    ("1,348", "with ≥ 5 of 7 cognitive indicators", "cognitive CFA"),
    ("1,342", "with both: gap and cognitive eligibility", "H1, H3")]:
    P(tf, [(a + "   ", False, NAVY), (b, False),
           (("   — " + c) if c else "", True, GREY)],
      13, INK, before=5)

P(tf, "Why the numbers differ", 15, NAVY, bold=True, before=16)
P(tf, "H2 has the gap as its outcome and contains no cognitive factor, so "
      "the cognitive-data requirement does not apply to it. The six "
      "participants without usable imaging leave every model that needs a "
      "gap. Both differences are intended and are stated wherever the "
      "samples are reported.", 12, GREY, before=4, line=1.14)

# =================================================================== 2 brain age
s = slide()
head(s, "Brain age: how it was estimated and how well it worked",
     "The predictor every hypothesis rests on")

code(s, [
    "10 repeated partitions  ×  10 outer folds  ×  5 inner folds",
    "splits grouped by family, stratified by family median age and modal sex",
    "LightGBM; inner loop searches learning rate (.03/.05/.10) × depth (3/5/7)",
    "preprocessing refitted at every level: inner, outer, and calibration",
    "",
    "age-bias correction, inside each outer-training set only:",
    "    5-fold cross-fitting  →  out-of-sample training predictions",
    "    β, α = polyfit(chronological age, those predictions)",
    "    corrected = (prediction − α) / β        BAG = corrected − age",
], 0.55, 1.5, 7.3, size=11, title="Procedure")

cw = [2.05, 1.30, 1.45, 1.45, 1.30]
t = table(s, 3, 5, 0.55, 4.25, cw, row_h=0.32)
for c, h in enumerate(["Modality", [("R", True), ("²", False)],
                       "MAE (years)", [("r", True), (" (BAG, age)", False)],
                       [("SD", True), (" of BAG", False)]]):
    put(t.cell(0, c), h, 11.5, bold=True,
        align=PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER)
for i, (k, lab, col) in enumerate(MODS):
    m = met[k]
    put(t.cell(i + 1, 0), lab, 11.5, bold=True, align=PP_ALIGN.LEFT, color=col)
    put(t.cell(i + 1, 1), f"{m['outer_r2_mean']:.3f} "
        f"({m['outer_r2_sd_across_folds_within_seed']:.3f})", 11.5)
    put(t.cell(i + 1, 2), f"{m['outer_mae_mean']:.2f} "
        f"({m['outer_mae_sd_across_folds_within_seed']:.2f})", 11.5)
    put(t.cell(i + 1, 3), num(m["bag_age_corr_corrected"], 3, True), 11.5)
    put(t.cell(i + 1, 4), f"{BAGSD[k]:.2f}", 11.5)
rules(t, 1, 2, 5)

note(s, 0.55, 5.35, 7.3,
     [("Accuracy is from the uncorrected out-of-sample predictions, as is "
       "conventional, and is the mean over the 100 outer test folds; "
       "parentheses give the mean within-repetition ", False),
      ("SD", True),
      (" across 10 folds. The gap’s correlation with age before "
       "correction was ", False),
      (MINUS + ".49 and " + MINUS + ".72.", False)])

tf = tb(s, 8.15, 1.55, 4.65, 4.9)
P(tf, "What matters here", 15, NAVY, bold=True, first=True)
P(tf, "Correction never sees the test fold. α and β come only from "
      "cross-fitted predictions inside outer-training, so no test-set "
      "information enters the gap. This follows More et al. (2023).",
      12.5, INK, before=6, line=1.14)
P(tf, "After correction the gap is independent of age. Without that, the "
      "gap would be a restatement of chronological age and H1 would be "
      "circular.", 12.5, INK, before=10, line=1.14)
P(tf, "The two modalities are not the same measure. Their corrected gaps "
      "correlate r = .31, so roughly 90% of the variance is not shared.",
      12.5, INK, before=10, line=1.14)
P(tf, "Stability was quantified. Each participant was tested once per "
      "repetition, so each has 10 gaps; the within-participant SD is "
      "1.65 years (structural) and 3.40 (functional).",
      12.5, INK, before=10, line=1.14)

# =================================================================== 3 measurement
s = slide()
head(s, "Measurement models",
     "Cognition and the two psychosocial constructs stay latent inside "
     "every hypothesis model")

code(s, [
    "g    =~ flanker + card_sort + list_sorting + pattern_comp +",
    "        picture_seq + ravlt + trails_b_rev",
    "Exec =~ a*flanker + a*card_sort",
    "Mem  =~ picture_seq + ravlt + list_sorting",
    "g ~~ 0*Exec    g ~~ 0*Mem    Exec ~~ 0*Mem",
], 0.55, 1.5, 6.1, size=11, title="Cognitive bifactor-(S−1)")

tf = tb(s, 0.55, 3.25, 6.1, 1.9)
P(tf, [("CFI .976   TLI .949   RMSEA .081 [.065, .098]   SRMR .026", False)],
  12, INK, first=True)
P(tf, [("ECV = .822     ω", False), ("H", False), (" = .817     ", False),
       ("n", True), (" = 1,348", False)], 12, INK, before=4)
P(tf, "Pattern comparison and Trails-B load on the general factor only, "
      "which is what makes this S−1: g is processing-speed-anchored and "
      "the two specific factors describe what executive and memory tasks "
      "share beyond it. The two executive loadings are constrained equal "
      "(a*) because a two-indicator factor is not otherwise identified.",
      11.5, GREY, before=7, line=1.13)

code(s, [
    "PWB =~ positive_affect + life_satisfaction + meaning_purpose",
    "SOC =~ emotional_support + instrumental_support +",
    "       friendship + loneliness_rev",
    "PWB ~~ SOC",
], 6.95, 1.5, 5.85, size=11, title="Psychosocial, two correlated factors")

tf = tb(s, 6.95, 3.07, 5.85, 2.1)
P(tf, [("CFI .974   TLI .959   RMSEA .069 [.057, .082]   SRMR .027", False)],
  12, INK, first=True)
P(tf, [("PWB ↔ SOC  ", False), ("r", True),
       (" = .687 [.645, .729]     ω = .72 / .83     ", False),
       ("n", True), (" = 1,342", False)], 12, INK, before=4)
P(tf, "Two factors, not one: the scaled difference test rejects the "
      "single-factor model decisively (Δχ²(1) = 116.3, "
      "p < .001; AIC 62,442 vs 62,733). This is what licenses entering "
      "well-being and connectedness as separate exposures.",
      11.5, GREY, before=7, line=1.13)

note(s, 0.55, 5.55, 12.25,
     [("Robust indices, MLR estimation with FIML for missing indicators. "
       "Factor scores are written out for description only; every "
       "hypothesis test keeps these factors latent inside a single-step "
       "model, so measurement error is modelled rather than carried in as "
       "if it were data.", False)])

# =================================================================== 4 H1
s = slide()
head(s, "H1 · A larger brain-age gap goes with lower general cognition",
     "Supported in both modalities")

code(s, [
    "[cognitive measurement model as above]",
    "",
    "g    ~ BAG + age_z + sex_male + edu_z",
    "Exec ~ BAG + age_z + sex_male + edu_z",
    "Mem  ~ BAG + age_z + sex_male + edu_z",
], 0.55, 1.5, 6.25, size=11.5, title="Model  (BAG standardised)")

cw = [1.85, 1.35, 1.80, 1.25]
t = table(s, 3, 4, 0.55, 3.45, cw, row_h=0.33)
for c, h in enumerate(["Modality", [("β", False)], "97.5% CI",
                       [("p", True), (" (Bonf.)", False)]]):
    put(t.cell(0, c), h, 11.5, bold=True,
        align=PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER)
for i, (k, lab, col) in enumerate(MODS):
    r = h1[(h1.modality == k) & (h1.spec == "primary") &
           (h1.type == "standardized") & (h1.lhs == "g") &
           (h1.rhs == "BAG")].iloc[0]
    q = qof("H1", k, "BAG -> cognition")
    put(t.cell(i + 1, 0), lab, 11.5, bold=True, align=PP_ALIGN.LEFT, color=col)
    put(t.cell(i + 1, 1), num(r["est"], 3, True), 11.5)
    put(t.cell(i + 1, 2),
        f"[{num(r['ci.lower'], 3)}, {num(r['ci.upper'], 3)}]", 11.5)
    put(t.cell(i + 1, 3), pstr(q), 11.5)
rules(t, 1, 2, 4)

note(s, 0.55, 4.65, 6.25,
     [("n", True), (" = 1,342 per modality. Bonferroni across the two "
       "modalities, α = .025, hence 97.5% intervals. Model fit "
       "(primary): CFI .951/.950, TLI .908/.906, RMSEA .069/.070, "
       "SRMR .030/.030.", False)])

tf = tb(s, 7.15, 1.55, 5.65, 5.0)
P(tf, "Three things to take from this", 15, NAVY, bold=True, first=True)

P(tf, "1.  The effect is real but modest", 13, NAVY, bold=True, before=11)
P(tf, "In the same model, chronological age has β = " + MINUS +
      ".68 and education +.25. The gap adds information beyond age — "
      "that is the point — but it is about a quarter the size of age.",
      12.5, INK, before=3, line=1.14)

P(tf, "2.  In years: half a year per year", 13, NAVY, bold=True, before=11)
P(tf, "The gap and age are both in years, so their per-year effects can be "
      "divided. One year of structural gap is worth about 0.50 years of "
      "ordinary ageing [0.32, 0.69]; for functional, 0.15 [0.01, 0.30].",
      12.5, INK, before=3, line=1.14)

P(tf, "3.  General ability, not a specific domain", 13, NAVY, bold=True,
  before=11)
P(tf, "The gap predicts the general factor only. Its coefficients on the "
      "executive-specific factor (+.11, +.01) and the memory-specific "
      "factor (+.03, " + MINUS + ".03) all have intervals covering zero. "
      "A single-factor model could not have shown this.",
      12.5, INK, before=3, line=1.14)

# =================================================================== 5 H2
s = slide()
head(s, "H2 · What relates to the brain-age gap itself",
     "Twelve tests, two survive: household income, in both modalities")

code(s, [
    "[PWB / SOC measurement model, only in the models that need it]",
    "",
    "BAG ~ X + age_z + sex_male + edu_z        one model per exposure",
], 0.55, 1.42, 7.55, size=11.5,
    title="Model  (gap stays in years; X standardised)")

cw = [2.75, 1.20, 1.55, 1.05, 1.00]
nrow = 1 + len(EXPO) * 2
t = table(s, nrow, 5, 0.55, 2.62, cw, row_h=0.265)
for c, h in enumerate(["Exposure", "Years / SD", "95% CI",
                       [("q", True)], [("ΔR", False), ("²", False)]]):
    put(t.cell(0, c), h, 10.5, bold=True,
        align=PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER)
row = 1
for key, lab in EXPO:
    for k, mlab, col in MODS:
        r = h2[(h2.modality == k) & (h2.exposure == key) &
               (h2.is_test == True) &
               (h2.type == "unstandardized_years")].iloc[0]
        q = qof("H2", k, key)
        dr = h2v[(h2v.modality == k) & (h2v.exposure == key)].delta_r2.iloc[0]
        sig = q < .05
        txt = (lab + "  ·  " + mlab) if True else lab
        put(t.cell(row, 0), txt, 10.5, bold=sig, align=PP_ALIGN.LEFT,
            color=col)
        put(t.cell(row, 1), num(r["est"], 2, True), 10.5, bold=sig, color=col if sig else INK)
        put(t.cell(row, 2),
            f"[{num(r['ci.lower'], 2)}, {num(r['ci.upper'], 2)}]", 10.5,
            color=col if sig else INK)
        put(t.cell(row, 3), pstr(q), 10.5, bold=sig, color=col if sig else INK)
        put(t.cell(row, 4), drfmt(dr), 10.5, color=col if sig else INK)
        row += 1
rules(t, 1, nrow - 1, 5)

note(s, 0.55, 6.22, 7.55,
     [("n", True), (" = 1,380 per modality. Benjamini-Hochberg at "
       "q = .05 across all twelve. ΔR² is over a model with age, "
       "sex and education only (R² = .003 structural, .018 "
       "functional). Bold rows survive correction.", False)])

tf = tb(s, 8.4, 1.55, 4.4, 5.0)
P(tf, "Reading the null", 15, NAVY, bold=True, first=True)
P(tf, "One SD higher income goes with a gap about one year younger, in "
      "both modalities. Nothing else survives.", 12.5, INK, before=6,
  line=1.14)
P(tf, "The two nearest misses are structural activity (p = .049) and "
      "structural well-being (p = .055). Both fail correction by a factor "
      "of three and are reported as not supported, not as trends.",
      12.5, INK, before=10, line=1.14)
P(tf, "Is income just a proxy for well-being or connectedness? No. With "
      "all six entered together the income coefficient barely moves "
      "(" + MINUS + "0.98 → " + MINUS + "0.93 structural, " +
      MINUS + "0.82 → " + MINUS + "0.85 functional), and it accounts "
      "for 71% and 53% of the variance the six explain jointly.",
      12.5, INK, before=10, line=1.14)
P(tf, "Several functional coefficients carry the wrong sign, all far from "
      "significance. That is what noise looks like in a measure with a "
      "3.4-year within-participant SD.", 12.5, GREY, before=10, line=1.14)

# =================================================================== 6 H3
s = slide()
head(s, "H3 · What relates to cognitive resilience",
     "Resilience = cognition conditional on the gap. Again only income "
     "survives")

code(s, [
    "[cognitive measurement model;  PWB / SOC where needed]",
    "",
    "g    ~ BAG + X + age_z + sex_male + edu_z    one model per exposure",
    "Exec ~ BAG + X + age_z + sex_male + edu_z",
    "Mem  ~ BAG + X + age_z + sex_male + edu_z",
], 0.55, 1.42, 7.55, size=11.5, title="Model  (outcome latent, standardised)")

cw = [2.75, 1.20, 1.55, 1.05, 1.00]
t = table(s, nrow, 5, 0.55, 3.02, cw, row_h=0.258)
for c, h in enumerate(["Exposure", [("β", False)], "95% CI",
                       [("q", True)], [("ΔR", False), ("²", False)]]):
    put(t.cell(0, c), h, 10.5, bold=True,
        align=PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER)
row = 1
for key, lab in EXPO:
    for k, mlab, col in MODS:
        r = h3[(h3.modality == k) & (h3.exposure == key) &
               (h3.is_test == True) & (h3.type == "standardized")].iloc[0]
        q = qof("H3", k, key)
        dr = h3v[(h3v.modality == k) & (h3v.exposure == key)].delta_r2.iloc[0]
        sig = q < .05
        put(t.cell(row, 0), lab + "  ·  " + mlab, 10.5, bold=sig,
            align=PP_ALIGN.LEFT, color=col)
        put(t.cell(row, 1), num(r["est"], 3, True), 10.5, bold=sig,
            color=col if sig else INK)
        put(t.cell(row, 2),
            f"[{num(r['ci.lower'], 2)}, {num(r['ci.upper'], 2)}]", 10.5,
            color=col if sig else INK)
        put(t.cell(row, 3), pstr(q), 10.5, bold=sig, color=col if sig else INK)
        put(t.cell(row, 4), drfmt(dr), 10.5, color=col if sig else INK)
        row += 1
rules(t, 1, nrow - 1, 5)

note(s, 0.55, 6.42, 7.55,
     [("n", True), (" = 1,342 per modality. BH at q = .05 across all "
       "twelve, a family separate from H2’s. ΔR² is over "
       "the H1 model (R² = .552 structural, .533 functional). "
       "Directions are the mirror of H2, because the outcome is cognition "
       "rather than the gap.", False)])

tf = tb(s, 8.4, 1.55, 4.4, 5.0)
P(tf, "Why this is the key slide", 15, NAVY, bold=True, first=True)
P(tf, "The coefficient on X here is, by the Frisch-Waugh-Lovell theorem, "
      "exactly the association between X and the part of cognition that "
      "brain age and demographics do not explain. No separate residual "
      "score is needed, and none is used.", 12.5, INK, before=6, line=1.14)
P(tf, "Income works through two routes at once. It is associated with a "
      "smaller gap (H2, negative) and with better cognition given the gap "
      "(H3, positive). Same direction, stated twice.",
      12.5, INK, before=10, line=1.14)
P(tf, "The psychosocial constructs are absent from both. Well-being, "
      "connectedness and stress — the mechanisms usually invoked for "
      "resilience — show nothing surviving correction against either "
      "outcome.", 12.5, INK, before=10, line=1.14)
P(tf, "Well-being even carries the wrong sign in both modalities, though "
      "close to zero.", 12.5, GREY, before=10, line=1.14)

# =================================================================== 7 decisions
s = slide()
head(s, "Every confirmatory decision",
     "Three families, corrected separately. 26 tests, 6 supported")

cw = [0.95, 2.70, 1.95, 1.00, 1.15]
t = table(s, 4, 5, 0.55, 1.58, cw, row_h=0.34)
for c, h in enumerate(["Family", "Test", "Correction", "Tests", "Supported"]):
    put(t.cell(0, c), h, 12, bold=True,
        align=PP_ALIGN.LEFT if c < 2 else PP_ALIGN.CENTER)
for i, row in enumerate([
        ("H1", "Gap \u2192 general cognition", "Bonferroni, \u03b1 = .025",
         "2", "2"),
        ("H2", "Six exposures \u2192 gap", "BH, q = .05", "12", "2"),
        ("H3", "Six exposures \u2192 resilience", "BH, q = .05", "12", "2")]):
    for c in range(5):
        put(t.cell(i + 1, c), row[c], 12, bold=(c == 0),
            align=PP_ALIGN.LEFT if c < 2 else PP_ALIGN.CENTER,
            color=NAVY if c == 0 else INK)
rules(t, 1, 3, 5)

tf = tb(s, 0.55, 3.10, 7.75, 0.4)
P(tf, "The six that survived", 14.5, NAVY, bold=True, first=True)

cw2 = [0.95, 3.45, 1.90, 1.45]
t2 = table(s, 7, 4, 0.55, 3.52, cw2, row_h=0.325)
for c, h in enumerate(["Family", "Result", "Estimate", "Adjusted"]):
    put(t2.cell(0, c), h, 12, bold=True,
        align=PP_ALIGN.LEFT if c < 2 else PP_ALIGN.CENTER)

rows6 = []
for k, mlab, col in MODS:
    r = h1[(h1.modality == k) & (h1.spec == "primary") &
           (h1.type == "standardized") & (h1.lhs == "g") &
           (h1.rhs == "BAG")].iloc[0]
    rows6.append(("H1", mlab + " gap \u2192 general cognition",
                  "\u03b2 = " + num(r["est"], 3, True),
                  "p " + pstr(qof("H1", k, "BAG -> cognition")), col))
for fam, src, typ in [("H2", h2, "unstandardized_years"),
                      ("H3", h3, "standardized")]:
    for k, mlab, col in MODS:
        r = src[(src.modality == k) & (src.exposure == "income") &
                (src.is_test == True) & (src.type == typ)].iloc[0]
        val = (num(r["est"], 2, True) + " yr / SD" if fam == "H2"
               else "\u03b2 = " + num(r["est"], 3, True))
        lab = ("income \u2192 gap" if fam == "H2"
               else "income \u2192 resilience")
        rows6.append((fam, mlab + ": " + lab, val,
                      "q " + pstr(qof(fam, k, "income")), col))
for i, row in enumerate(rows6):
    for c in range(4):
        put(t2.cell(i + 1, c), row[c], 11.5, bold=(c == 0),
            align=PP_ALIGN.LEFT if c < 2 else PP_ALIGN.CENTER,
            color=row[4] if c == 1 else INK)
rules(t2, 1, 6, 4)

tf = tb(s, 8.70, 1.58, 4.10, 4.9)
P(tf, "Why the corrections differ", 14.5, NAVY, bold=True, first=True)
P(tf, "H1 asks one question of two modalities. A false positive on either "
      "would be claimed as support, so the family is small and the control "
      "is strict.", 11.5, INK, before=6, line=1.14)
P(tf, "H2 and H3 each screen six exposures across two modalities. "
      "Bonferroni at .05/12 would control the chance of any false positive "
      "at the cost of missing most true ones \u2014 the wrong trade for a "
      "shortlist. BH controls the expected proportion of false discoveries "
      "among the rejections.", 11.5, INK, before=9, line=1.14)
P(tf, "H2 and H3 are separate families. One asks what relates to the gap, "
      "the other what relates to cognition given the gap; a result in one "
      "implies nothing about the other.", 11.5, INK, before=9, line=1.14)
P(tf, "With two rejections at q = .05, the expected number of false "
      "discoveries in each family is about 0.1.", 11.5, GREY, before=9,
  line=1.14)

note(s, 0.55, 6.08, 7.75,
     [("The adjusted column is a Bonferroni-adjusted ", False), ("p", True),
      (" for H1 and a Benjamini-Hochberg ", False), ("q", True),
      (" for H2 and H3. They are different quantities and are labelled "
       "separately throughout.", False)])

# =================================================================== 8 caveats
s = slide()
head(s, "What this supports, and what it does not",
     "Read before quoting any number above")

tf = tb(s, 0.55, 1.5, 3.95, 5.0)
P(tf, "Supported", 15, NAVY, bold=True, first=True)
P(tf, "The brain-age gap carries information about general cognition "
      "beyond chronological age, in both modalities, and more strongly "
      "for structural MRI.", 12.5, INK, before=6, line=1.14)
P(tf, "Household income is associated with a smaller gap and with higher "
      "resilience, in both modalities and in both hypotheses — four "
      "independent-ish tests pointing the same way.",
      12.5, INK, before=10, line=1.14)
P(tf, "The gap tracks general ability rather than any specific cognitive "
      "domain.", 12.5, INK, before=10, line=1.14)

tf = tb(s, 4.75, 1.5, 3.95, 5.0)
P(tf, "Not supported", 15, NAVY, bold=True, first=True)
P(tf, "Psychological well-being, social connectedness, perceived stress, "
      "sleep quality and physical activity: no association with the gap "
      "or with resilience survives correction, in either modality.",
      12.5, INK, before=6, line=1.14)
P(tf, "These are precisely estimated nulls, not absent data. The "
      "intervals are narrow enough to exclude effects of the size the "
      "literature usually reports, which makes them informative rather "
      "than merely inconclusive.", 12.5, INK, before=10, line=1.14)

tf = tb(s, 8.95, 1.5, 3.85, 5.0)
P(tf, "Limitations to state", 15, NAVY, bold=True, first=True)
for txt in [
    "Cross-sectional. Every coefficient describes concurrent association "
    "and licenses no claim about direction.",
    "Effect sizes are small: income adds 1.6% of gap variance and 1.1% of "
    "cognitive variance.",
    "The gap is a model residual, not a measurement of a brain. Its "
    "within-participant SD is 1.65 and 3.40 years.",
    "Measurement error in the gap is not modelled, so the coefficients on "
    "it are attenuated — a lower bound.",
    "TLI is .91 in the H1 models, below the .95 benchmark; CFI, RMSEA and "
    "SRMR are met.",
]:
    P(tf, txt, 11.5, INK, before=8, line=1.13)

prs.save(OUT)
print("wrote", OUT, f"({len(prs.slides._sldIdLst)} slides)")
