"""The five figures as one deck, in the order that teaches the study.

Each slide is one figure and one question. The question is the title, the
answer is the subtitle, and the figure is the evidence. A reader who sees
only the titles should still come away with the argument:

    Can brain age be estimated?                yes, structural better
    How much do people differ?                 substantially
    Does the gap matter for cognition?         yes, both modalities
    What is associated with gap and resilience? income, and only income
    What does the study amount to?             the graphical abstract

Figures are built with FIG_NO_TITLE so they carry no heading of their
own; the slide supplies it in the deck's typography.
"""
import os
import sys

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
from PIL import Image

FIGD = sys.argv[1]
OUT = sys.argv[2]

NAVY = RGBColor(0x16, 0x32, 0x4F)
GREY = RGBColor(0x5A, 0x65, 0x70)
INK = RGBColor(0x1A, 0x1A, 0x1A)
BLUE = RGBColor(0x15, 0x49, 0x7B)
PURPLE = RGBColor(0x7B, 0x2D, 0x8E)

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]


def tf_at(s, x, y, w, h):
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


def picture(s, name, left, top, max_w, max_h):
    path = os.path.join(FIGD, name)
    if not os.path.exists(path):
        raise SystemExit("missing figure: " + path)
    iw, ih = Image.open(path).size
    k = min(max_w / iw, max_h / ih)
    w, h = iw * k, ih * k
    s.shapes.add_picture(path, Inches(left + (max_w - w) / 2),
                         Inches(top + (max_h - h) / 2), Inches(w), Inches(h))


def fig_slide(number, question, answer, png, footnote, fig_top=1.62,
              fig_h=5.05):
    s = prs.slides.add_slide(BLANK)
    tf = tf_at(s, 0.6, 0.38, 12.1, 1.1)
    P(tf, [(f"Figure {number}   ", False, GREY), (question, False, NAVY)],
      25, NAVY, font="Cambria", bold=True, first=True)
    P(tf, answer, 14.5, GREY, before=5)
    picture(s, png, 0.5, fig_top, 12.33, fig_h)
    if footnote:
        tf = tf_at(s, 0.6, 6.82, 12.1, 0.5)
        P(tf, footnote, 11.5, GREY, first=True, line=1.12)
    return s


# ------------------------------------------------------------------ title
s = prs.slides.add_slide(BLANK)
tf = tf_at(s, 1.1, 2.15, 11.1, 2.6)
P(tf, "Same brain age, different cognition", 40, NAVY, font="Cambria",
  bold=True, first=True)
P(tf, "Operationalising cognitive resilience as deviation from brain age",
  20, GREY, before=10)
P(tf, [("AABC / HCP-Aging Release 2   ·   ", False),
       ("n", True), (" = 1,380 adults aged 36–90", False)],
  15, GREY, before=18)

tf = tf_at(s, 1.1, 4.62, 11.1, 2.55)
P(tf, "Each figure answers one question, in order", 13.5, NAVY,
  bold=True, first=True)
for i, q in enumerate([
        "Can brain age be estimated, and what did the correction do?",
        "How much do people actually differ in the gap?",
        "Does the gap matter for cognition?",
        "Two people, the same brain age, different cognition",
        "How is cognitive resilience computed?",
        "What is associated with the brain-age gap?",
        "What is associated with cognitive resilience?",
        "What does the whole study amount to?"]):
    P(tf, [(f"{i + 1}.  ", False, NAVY), (q, False)], 12.5, INK, before=2)

# ------------------------------------------------------------------ 1
fig_slide(
    1, "Can brain age be estimated from MRI?",
    "Yes. Structural MRI predicted chronological age more accurately than "
    "functional connectivity, and the bias correction largely removed the "
    "gap’s dependence on age.",
    "fig1_brainage_presentation.png",
    [("Out-of-sample predictions from strictly nested, family-grouped "
      "cross-validation: 10 repeated partitions × 10 outer × 5 "
      "inner folds. Accuracy is the mean over the 100 outer test folds, "
      "computed before correction. Correction coefficients were estimated "
      "only within outer-training data.", False)],
    fig_top=1.52, fig_h=5.20)

# ------------------------------------------------------------------ 5 -> 2
fig_slide(
    2, "How much do people differ in the gap?",
    "Substantially. The middle half of participants span about ±5 "
    "years (structural), and the tails reach ±30.",
    "fig9_bag_violin_presentation.png",
    [("Every participant is drawn. Between-person SD is the spread across "
      "people; it is not the uncertainty of one person’s value, which "
      "is 1.65 years (structural) and 3.40 (functional) across the ten "
      "repeated partitions. Open rings mark the two participants beyond "
      "four SDs; none was excluded.", False)],
    fig_top=1.55, fig_h=5.10)

# ------------------------------------------------------------------ 2 -> 3
fig_slide(
    3, "Does the gap matter for cognition?",
    "Yes, in both modalities. A larger gap went with lower general "
    "cognitive performance, over and above chronological age.",
    "fig2_h1_presentation.png",
    [("Adjusted for chronological age, sex and education. Intervals are "
      "97.5%, reflecting the Bonferroni correction across the two "
      "modalities. ", False), ("n", True),
     (" = 1,342 per modality. The two coefficients were not tested against "
      "each other, so the figure does not claim the structural "
      "association is larger.", False)],
    fig_top=1.95, fig_h=4.35)

# ------------------------------------- resilience, the core idea, in two
fig_slide(
    4, "Two people, the same brain age, different cognition",
    "The brain-age gap does not fix cognitive performance. What is left "
    "over is what this thesis measures.",
    "fig5a_resilience_concept_presentation.png",
    [("Every participant is drawn. The shaded band holds the 148 people "
      "whose structural gap falls between +6 and +10 years — the same "
      "brain age, in any practical sense — and their cognition spans "
      "4.5 standard deviations. The two ringed participants are the "
      "extremes of that band.", False)],
    fig_top=1.62, fig_h=4.85)

fig_slide(
    5, "How cognitive resilience is computed",
    "Predict each person’s cognition from their brain and background, "
    "then keep the difference between that prediction and what they "
    "actually scored.",
    "fig5b_resilience_recipe_presentation.png",
    [("The two people shown are real participants from the band in the "
      "previous figure, not an illustration. Because resilience is defined "
      "as the part of cognition the gap does not explain, it is "
      "uncorrelated with the gap by construction (r = .000) — that is "
      "a property of the definition, not a result. A positive score means "
      "the person scored better than the model predicted; it does not mean "
      "they are protected.", False)],
    fig_top=1.62, fig_h=4.80)

# ------------------------------------------------------------------ 3 -> 4
fig_slide(
    6, "What is associated with the brain-age gap?",
    "Of six lifestyle and psychosocial exposures, only household income "
    "survived correction \u2014 in both modalities.",
    "fig6a_h2_gap_presentation.png",
    [("Each exposure was tested in its own model, adjusted for age, sex "
      "and education. Benjamini-Hochberg at q = .05 across the twelve "
      "tests of this family. One open marker \u2014 structural physical "
      "activity \u2014 has an interval excluding zero but did not survive "
      "correction, which is why the decision is carried by the marker and "
      "not by the interval. ", False), ("n", True),
     (" = 1,380 per modality.", False)],
    fig_top=1.62, fig_h=4.85)

fig_slide(
    7, "What is associated with cognitive resilience?",
    "The same answer, from a different question: only household income, "
    "again in both modalities.",
    "fig6b_h3_resilience_presentation.png",
    [("Same six exposures, same one-model-per-exposure design, now "
      "adjusted for the brain-age gap as well. A separate "
      "Benjamini-Hochberg family of twelve, because this question is not "
      "the previous one. Functional physical activity is the open marker "
      "whose interval excludes zero. ", False), ("n", True),
     (" = 1,342 per modality. Coefficients here are standardised and are "
      "not comparable in size with the years on the previous figure.",
      False)],
    fig_top=1.45, fig_h=5.55)

# ------------------------------------------------------------------ 4 -> 5
fig_slide(
    8, "What does the study amount to?",
    "Four findings, in order.",
    "fig4_summary_presentation.png",
    [("All data are cross-sectional and collected at a single visit. Every "
      "relationship shown is an association between concurrently measured "
      "variables and licenses no claim about direction: income may support "
      "cognition, and cognition may support income.", False)],
    fig_top=1.62, fig_h=4.95)

prs.save(OUT)
print("wrote", OUT, f"({len(prs.slides._sldIdLst)} slides)")
