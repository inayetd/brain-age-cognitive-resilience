"""Results slides: two APA 7 tables and the seven brain-age figures.

Typography follows the other decks in presentations/ (Cambria headings in
#16324F, Calibri body, #5A6570 for secondary text) so this can be dropped in.

APA 7 table rules applied: horizontal rules only, above the header, under the
column heads and under the body; a short rule under each spanner head; no
vertical lines; no leading zero on bounded statistics; a real minus sign; the
note introduced by an italic "Note.".

The figure images used here are the presentation style with the baked-in title
suppressed, because the slide supplies the heading instead.
"""
import os
import sys

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt, Emu
from pptx.oxml.ns import qn
from lxml import etree
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

MINUS = '−'


def slide():
    return prs.slides.add_slide(BLANK)


def textbox(s, x, y, w, h):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    return tf


def put(tf, text, size, color, font='Calibri', bold=False, italic=False,
        align=PP_ALIGN.LEFT, first=False, space_before=0):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    p.space_before = Pt(space_before)
    r = p.add_run()
    r.text = text
    r.font.size = Pt(size)
    r.font.color.rgb = color
    r.font.name = font
    r.font.bold = bold
    r.font.italic = italic
    return p


def rich(tf, parts, size, color, font='Calibri', align=PP_ALIGN.LEFT,
         first=False, space_before=0):
    """parts: list of (text, italic) or (text, italic, 'sup')."""
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    p.space_before = Pt(space_before)
    for part in parts:
        t, it = part[0], part[1]
        sup = len(part) > 2 and part[2] == 'sup'
        r = p.add_run()
        r.text = t
        r.font.size = Pt(size)
        r.font.color.rgb = color
        r.font.name = font
        r.font.italic = it
        if sup:
            r.font._rPr.set('baseline', '30000')
            r.font.size = Pt(size * 0.8)
    return p


def heading(s, title, subtitle=None):
    tf = textbox(s, 0.62, 0.42, 12.1, 1.1)
    put(tf, title, 26, NAVY, font='Cambria', bold=True, first=True)
    if subtitle:
        put(tf, subtitle, 13.5, GREY, space_before=5)


# ------------------------------------------------------------------ borders
def set_border(cell, edge, on):
    """edge in {'L','R','T','B'}. APA tables carry no vertical rules, so the
    L and R edges are always written as explicit noFill: leaving them unset
    lets the inherited table style draw a grid."""
    tc = cell._tc
    tcPr = tc.get_or_add_tcPr()
    tag = qn('a:ln' + edge)
    for old in tcPr.findall(tag):
        tcPr.remove(old)
    ln = etree.SubElement(tcPr, tag)
    ln.set('w', '12700')          # 1 pt
    ln.set('cap', 'flat')
    ln.set('cmpd', 'sng')
    ln.set('algn', 'ctr')
    if on:
        fill = etree.SubElement(ln, qn('a:solidFill'))
        clr = etree.SubElement(fill, qn('a:srgbClr'))
        clr.set('val', '000000')
    else:
        etree.SubElement(ln, qn('a:noFill'))
    # order matters to the schema: lnL, lnR, lnT, lnB must follow the fill
    for e in ('L', 'R', 'T', 'B'):
        t = qn('a:ln' + e)
        found = tcPr.findall(t)
        for f in found:
            tcPr.remove(f)
            tcPr.append(f)


def strip_fill(cell):
    tcPr = cell._tc.get_or_add_tcPr()
    for t in ('a:solidFill', 'a:noFill', 'a:gradFill', 'a:blipFill'):
        for old in tcPr.findall(qn(t)):
            tcPr.remove(old)
    etree.SubElement(tcPr, qn('a:noFill'))


def make_table(s, rows, cols, x, y, w, col_w, row_h):
    shp = s.shapes.add_table(rows, cols, Inches(x), Inches(y),
                             Inches(w), Inches(row_h * rows))
    tbl = shp.table
    # Kill the default banded theme: an APA table has no shading.
    tblPr = tbl._tbl.find(qn('a:tblPr'))
    tblPr.set('firstRow', '0')
    tblPr.set('bandRow', '0')
    for st in tblPr.findall(qn('a:tableStyleId')):
        tblPr.remove(st)
    for i, cw in enumerate(col_w):
        tbl.columns[i].width = Inches(cw)
    for r in range(rows):
        tbl.rows[r].height = Inches(row_h)
    for r in range(rows):
        for c in range(cols):
            cell = tbl.cell(r, c)
            strip_fill(cell)
            for e in ('L', 'R', 'T', 'B'):
                set_border(cell, e, False)
            cell.margin_left = Inches(0.06)
            cell.margin_right = Inches(0.06)
            cell.margin_top = Inches(0.03)
            cell.margin_bottom = Inches(0.03)
            cell.vertical_anchor = MSO_ANCHOR.BOTTOM
    return tbl


def fill(cell, parts, size=15, bold=False, align=PP_ALIGN.CENTER, color=INK):
    tf = cell.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = align
    if isinstance(parts, str):
        parts = [(parts, False)]
    for part in parts:
        t, it = part[0], part[1]
        sup = len(part) > 2 and part[2] == 'sup'
        r = p.add_run()
        r.text = t
        r.font.size = Pt(size * 0.8 if sup else size)
        r.font.name = 'Calibri'
        r.font.italic = it
        r.font.bold = bold
        r.font.color.rgb = color
        if sup:
            r.font._rPr.set('baseline', '30000')


def reorder(tbl):
    """CT_TableCellProperties fixes the order of its children: the four line
    elements first, the fill last. Written in any other order the properties
    are silently dropped, which is invisible until a renderer ignores the
    rules. Called once after every cell is finished."""
    order = ['a:lnL', 'a:lnR', 'a:lnT', 'a:lnB', 'a:lnTlToBr', 'a:lnBlToTr',
             'a:cell3D', 'a:noFill', 'a:solidFill', 'a:gradFill', 'a:blipFill',
             'a:pattFill', 'a:grpFill', 'a:headers', 'a:extLst']
    for row in tbl.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            kids = list(tcPr)
            for k in kids:
                tcPr.remove(k)
            for name in order:
                for k in kids:
                    if k.tag == qn(name):
                        tcPr.append(k)


def rule_between(tbl, r_above, cols):
    """Draw a rule on the shared edge between r_above and the row under it.

    Both sides have to be written. A cell's top edge and the cell above it
    share one line, and a renderer that finds lnB solid above and lnT noFill
    below resolves the conflict in favour of the lower cell, so setting only
    lnB produces nothing at all.
    """
    for c in cols:
        set_border(tbl.cell(r_above, c), 'B', True)
        if r_above + 1 < len(tbl.rows):
            set_border(tbl.cell(r_above + 1, c), 'T', True)


def no_verticals(tbl):
    """Re-assert the noFill verticals after any merge, which can reintroduce
    the inherited style's grid on the cells it rewrites."""
    for row in tbl.rows:
        for cell in row.cells:
            for e in ('L', 'R'):
                set_border(cell, e, False)


R2 = [('R', True), ('2', False, 'sup')]

# ==================================================================== slide 1
s = slide()
heading(s, 'Table 1',
        'What the age-prediction models achieve, and what the gap they '
        'produce looks like')

tf = textbox(s, 0.62, 1.62, 12.1, 0.4)
put(tf, 'Brain-Age Prediction Accuracy and Properties of the Brain-Age Gap',
    16, INK, italic=True, first=True)

# 8 columns across 11.9 in
cw1 = [1.75, 0.95, 1.55, 1.70, 1.35, 1.35, 1.65, 1.60]
t = make_table(s, 4, 8, 0.62, 2.18, sum(cw1), cw1, 0.46)

t.cell(0, 4).merge(t.cell(0, 5))
t.cell(0, 6).merge(t.cell(0, 7))
for c in (0, 1, 2, 3):
    t.cell(0, c).merge(t.cell(1, c))

fill(t.cell(0, 0), 'Modality', bold=True, align=PP_ALIGN.LEFT)
fill(t.cell(0, 1), [('n', True)], bold=True)
fill(t.cell(0, 2), R2, bold=True)
fill(t.cell(0, 3), 'MAE (years)', bold=True)
fill(t.cell(0, 4), [('r', True), (' (BAG, age)', False)], bold=True)
fill(t.cell(0, 6), [('SD', True), (' of BAG (years)', False)], bold=True)
fill(t.cell(1, 4), 'Before', bold=True)
fill(t.cell(1, 5), 'After', bold=True)
fill(t.cell(1, 6), 'Between persons', bold=True)
fill(t.cell(1, 7), 'Within person', bold=True)

body1 = [
    ('Structural', '1,380', '.82 (.03)', '5.35 (0.39)',
     MINUS + '.49', '.02', '7.42', '1.65', BLUE),
    ('Functional', '1,380', '.69 (.04)', '7.20 (0.39)',
     MINUS + '.72', '.08', '9.94', '3.40', PURPLE),
]
for i, row in enumerate(body1):
    r = 2 + i
    fill(t.cell(r, 0), row[0], bold=True, align=PP_ALIGN.LEFT, color=row[8])
    for c in range(1, 8):
        fill(t.cell(r, c), row[c])

no_verticals(t)
for c in range(8):
    set_border(t.cell(0, c), 'T', True)          # rule above the table
rule_between(t, 0, (4, 6))                       # short rules under the spanners
rule_between(t, 1, range(8))                     # rule under the column heads
for c in range(8):
    set_border(t.cell(3, c), 'B', True)          # rule under the body
reorder(t)

tf = textbox(s, 0.62, 4.30, 11.9, 2.4)
rich(tf, [('Note. ', True),
          ('Accuracy is computed from the ', False), ('uncorrected', True),
          (' out-of-sample predictions and is the mean across the 100 outer '
           'test folds (10 folds × 10 repetitions); parentheses give the '
           'mean within-repetition ', False), ('SD', True),
          (' across 10 folds. ', False), ('r', True),
          (' (BAG, age) is the correlation of the gap with chronological age '
           'before and after the cross-fitted bias correction. '
           'Between persons = spread of the gap across people; '
           'within person = stability of one person’s gap across the 10 '
           'repetitions. BAG = brain-age gap; MAE = mean absolute error.',
           False)],
     13, GREY, first=True)

# ==================================================================== slide 2
s = slide()
heading(s, 'Table 2',
        'Accuracy of the predicted age after correction — a calibration '
        'check, not the accuracy of the models')

tf = textbox(s, 0.62, 1.62, 12.1, 0.4)
put(tf, 'Calibration Diagnostic: Accuracy of the Bias-Corrected Predicted Age',
    16, INK, italic=True, first=True)

cw2 = [2.40, 1.60, 2.10, 2.10, 1.60]
t = make_table(s, 3, 5, 0.62, 2.18, sum(cw2), cw2, 0.50)

fill(t.cell(0, 0), 'Modality', bold=True, align=PP_ALIGN.LEFT)
fill(t.cell(0, 1), R2, bold=True)
fill(t.cell(0, 2), 'MAE (years)', bold=True)
fill(t.cell(0, 3), 'RMSE (years)', bold=True)
fill(t.cell(0, 4), [('b', True)], bold=True)

body2 = [('Structural', '.78', '5.80', '7.42', '1.01', BLUE),
         ('Functional', '.61', '7.90', '9.94', '1.05', PURPLE)]
for i, row in enumerate(body2):
    r = 1 + i
    fill(t.cell(r, 0), row[0], bold=True, align=PP_ALIGN.LEFT, color=row[5])
    for c in range(1, 5):
        fill(t.cell(r, c), row[c])

no_verticals(t)
for c in range(5):
    set_border(t.cell(0, c), 'T', True)          # rule above the table
    set_border(t.cell(2, c), 'B', True)          # rule under the body
rule_between(t, 0, range(5))                     # rule under the column heads
reorder(t)

tf = textbox(s, 0.62, 3.95, 11.9, 2.6)
rich(tf, [('Note. ', True),
          ('Values are pooled over all participants and all 10 repetitions, '
           'not averaged over folds, because the corrected predictions were '
           'never scored fold by fold. They are therefore ', False),
          ('not', True),
          (' comparable with Table 1 and do not represent the predictive '
           'accuracy of the models. ', False),
          ('b', True),
          (' is the slope of corrected predicted age on chronological age; '
           'b ≈ 1 means the age-related bias has been removed. '
           'MAE = mean absolute error; RMSE = root mean square error.',
           False)],
     13, GREY, first=True)

# ==================================================================== slide 3
# The PCA / no-PCA comparison. Two questions with different answers, so the
# table carries both: accuracy, which is a property of the model, and the gap,
# which is a property of a person.
s = slide()
heading(s, 'Table 3',
        'Does the dimension reduction cost anything, and do the two '
        'pipelines agree about individuals?')

tf = textbox(s, 0.62, 1.62, 12.1, 0.4)
put(tf, 'Functional Pipeline With and Without Principal Component Analysis',
    16, INK, italic=True, first=True)

cw3 = [2.20, 2.00, 1.45, 1.75, 1.45, 1.60, 1.45]
t = make_table(s, 4, 7, 0.62, 2.18, sum(cw3), cw3, 0.46)

t.cell(0, 5).merge(t.cell(0, 6))
for c in (0, 1, 2, 3, 4):
    t.cell(0, c).merge(t.cell(1, c))

fill(t.cell(0, 0), 'Pipeline', bold=True, align=PP_ALIGN.LEFT)
fill(t.cell(0, 1), 'Features in the model', bold=True)
fill(t.cell(0, 2), R2, bold=True)
fill(t.cell(0, 3), 'MAE (years)', bold=True)
fill(t.cell(0, 4), [('r', True), (' (BAG, age)', False)], bold=True)
fill(t.cell(0, 5), [('SD', True), (' of BAG (years)', False)], bold=True)
fill(t.cell(1, 5), 'Between persons', bold=True)
fill(t.cell(1, 6), 'Within person', bold=True)

body3 = [('Functional (PCA)', '100 components', '.69 (.03)', '7.20 (0.39)',
          '.08', '9.94', '3.40', PURPLE),
         ('Functional (no PCA)', '71,631 edges', '.69 (.03)', '7.17 (0.38)',
          '.05', '10.41', '2.63', RGBColor(0x11, 0x7A, 0x65))]
for i, row in enumerate(body3):
    r = 2 + i
    fill(t.cell(r, 0), row[0], bold=True, align=PP_ALIGN.LEFT, color=row[7])
    for c in range(1, 7):
        fill(t.cell(r, c), row[c])

no_verticals(t)
for c in range(7):
    set_border(t.cell(0, c), 'T', True)
    set_border(t.cell(3, c), 'B', True)
rule_between(t, 0, (5,))
rule_between(t, 1, range(7))
reorder(t)

tf = textbox(s, 0.62, 4.30, 11.9, 2.4)
rich(tf, [('Note. ', True),
          ('Both pipelines were run over the same 10 repeated partitions with '
           'the same folds, so the fold-level values are paired. Compared that '
           'way across the 10 repetitions, PCA loses 0.008 in ', False),
          ('R', True), ('2', False, 'sup'),
          (' (', False), ('t', True),
          ('(9) = −4.51, p = .001) and gains 0.04 years of MAE '
           '(', False), ('t', True), ('(9) = 1.67, p = .129): a difference '
           'that is detectable but of no practical size. The two corrected '
           'gaps nonetheless correlate only ', False), ('r', True),
          (' = .72, the median person differs by 5.0 years and half the '
           'sample by more than 5 years, so the choice matters for who is '
           'called brain-older even though it does not matter for accuracy. '
           'PCA is the pre-specified primary pipeline; no PCA is reported as '
           'a sensitivity analysis. BAG = brain-age gap; MAE = mean absolute '
           'error.', False)],
     13, GREY, first=True)

# ==================================================================== slide 4
# Range and tails of the gap. Percentiles are deliberately left out: they
# answer a question nobody asked out loud, and the minimum, maximum and SD
# already describe the spread.
s = slide()
heading(s, 'Table 4',
        'How far the gap spreads, and how many participants sit far enough '
        'out to matter')

tf = textbox(s, 0.62, 1.62, 12.1, 0.4)
put(tf, 'Distribution of the Corrected Brain-Age Gap', 16, INK, italic=True,
    first=True)

cw4 = [1.95, 0.90, 0.95, 0.95, 1.15, 1.15, 1.15, 1.00, 1.30, 1.40]
t = make_table(s, 4, 10, 0.62, 2.18, sum(cw4), cw4, 0.46)

heads4 = ['Modality', [('n', True)], [('M', True)], [('SD', True)],
          'Min', 'Max', 'Range', 'Skew', 'Kurtosis', [('n', True),
                                                      (' with |z| > 4', False)]]
for c, h in enumerate(heads4):
    fill(t.cell(0, c), h, bold=True,
         align=PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER)

body4 = [
    ('Structural', '1,380', '0.00', '7.42', '−21.99', '31.46', '53.45',
     '0.06', '0.56', '2', BLUE),
    ('Functional (PCA)', '1,380', '−0.14', '9.94', '−32.57', '34.40',
     '66.96', '0.16', '0.17', '0', PURPLE),
    ('Functional (no PCA)', '1,380', '−0.01', '10.41', '−33.84',
     '35.30', '69.14', '0.04', '−0.09', '0',
     RGBColor(0x11, 0x7A, 0x65)),
]
for i, row in enumerate(body4):
    r = 1 + i
    fill(t.cell(r, 0), row[0], bold=True, align=PP_ALIGN.LEFT, color=row[10])
    for c in range(1, 10):
        fill(t.cell(r, c), row[c])

no_verticals(t)
for c in range(10):
    set_border(t.cell(0, c), 'T', True)
    set_border(t.cell(3, c), 'B', True)
rule_between(t, 0, range(10))
reorder(t)

tf = textbox(s, 0.62, 4.15, 11.9, 2.6)
rich(tf, [('Note. ', True),
          ('All values in years except skew, kurtosis and the count. '
           'Kurtosis is excess kurtosis. ', False),
          ('z', True),
          (' is standardised within modality, so |', False), ('z', True),
          ('| > 4 means four standard deviations from that modality’s own '
           'mean; under normality the expected count at ', False),
          ('n', True),
          (' = 1,380 is 0.09 per modality. The two structural cases are both '
           'positive (+31.5 and +30.3 years, in participants aged 43 and 56) '
           'and sit at the end of a continuous tail rather than detached from '
           'it. No participant was excluded on the basis of these values. The '
           'wider range belongs to the functional models, as their lower '
           'accuracy implies; the only cases beyond |', False), ('z', True),
          ('| > 4 belong to the structural model, whose narrower and more '
           'peaked distribution makes the same absolute deviation larger in '
           'standardised terms. No participant is extreme in more than one '
           'modality.', False)],
     13, GREY, first=True)

# ============================================================ figure slides
FIGS = [
    ('fig1_prediction_accuracy',
     'Accuracy of age prediction',
     'Predicted against chronological age, uncorrected predictions'),
    ('fig7_prediction_corrected',
     'Accuracy of age prediction, before and after bias correction',
     'Top row uncorrected, bottom row corrected'),
    ('fig2_bias_correction',
     'Brain-age gap against chronological age',
     'Before and after the cross-fitted bias correction'),
    ('fig9_bag_violin',
     'Distribution of the brain-age gap',
     'Every participant shown; the gap is a position in this distribution, '
     'not a statement about one person'),
    ('fig3_gap_distribution',
     'Distribution of the brain-age gap',
     'Spread across people, and the gap’s stability within a person'),
    ('fig4_error_by_age_band',
     'Mean brain-age gap by age band',
     'Before and after correction'),
    # fig8 supersedes fig5 here: it carries the same agreement panel and adds
    # the paired accuracy comparison, so fig5 would only repeat it.
    ('fig8_pca_vs_nopca',
     'Functional pipeline with and without PCA',
     'Same accuracy, different people'),
    # fig6 (gap by site and sex) is deliberately not a slide. It is still
    # produced for the thesis, where it supports a limitation, but it is not
    # part of what gets presented.
]

for name, title, sub in FIGS:
    path = os.path.join(FIGD, name + '_presentation.png')
    if not os.path.exists(path):
        raise SystemExit('missing ' + path)
    s = slide()
    heading(s, title, sub)
    iw, ih = Image.open(path).size
    top, bottom = 1.72, 0.35
    avail_h = 7.5 - top - bottom
    avail_w = 12.1
    scale = min(avail_w / iw, avail_h / ih)
    w, h = iw * scale, ih * scale
    s.shapes.add_picture(path, Inches((13.333 - w) / 2), Inches(top),
                         Inches(w), Inches(h))

prs.save(OUT)
print('wrote ' + OUT + f'  ({len(prs.slides.__iter__.__self__._sldIdLst)} slides)')
