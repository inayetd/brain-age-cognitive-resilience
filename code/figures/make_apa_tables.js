// APA 7 tables for the brain-age results.
// Horizontal rules only: above the header, under the column heads, under the body.
// Spanner heads get a rule under the columns they span. No vertical lines anywhere.
const d = require('docx');
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  WidthType, BorderStyle, AlignmentType, VerticalAlign, PageOrientation,
} = d;

const FONT = 'Times New Roman';
const BODY = 24;   // half-points -> 12 pt
const TBL = 20;    // half-points -> 10 pt

const LINE = { style: BorderStyle.SINGLE, size: 6, color: '000000' };
const NONE = { style: BorderStyle.NONE, size: 0, color: 'FFFFFF' };

// Cell borders: pass which edges carry a rule.
function bord({ top = false, bottom = false } = {}) {
  return {
    top: top ? LINE : NONE,
    bottom: bottom ? LINE : NONE,
    left: NONE,
    right: NONE,
  };
}

function runs(spec) {
  // spec: string, or array of {t, i (italic), sup, sub}
  if (typeof spec === 'string') spec = [{ t: spec }];
  return spec.map(s => new TextRun({
    text: s.t,
    italics: !!s.i,
    superScript: !!s.sup,
    subScript: !!s.sub,
    font: FONT,
    size: s.size || TBL,
  }));
}

function cell(spec, { width, span = 1, rowSpan = 1, align = AlignmentType.CENTER,
                      top = false, bottom = false } = {}) {
  return new TableCell({
    width: { size: width, type: WidthType.DXA },
    columnSpan: span,
    rowSpan: rowSpan,
    borders: bord({ top, bottom }),
    verticalAlign: VerticalAlign.BOTTOM,
    margins: { top: 60, bottom: 60, left: 80, right: 80 },
    children: [new Paragraph({
      alignment: align,
      spacing: { before: 0, after: 0, line: 240 },
      children: runs(spec),
    })],
  });
}

function para(children, { spacing = { before: 0, after: 0, line: 240 }, align } = {}) {
  return new Paragraph({ alignment: align, spacing, children });
}

function tableNumber(n) {
  return para([new TextRun({ text: `Table ${n}`, bold: true, font: FONT, size: BODY })],
              { spacing: { before: 0, after: 0, line: 480 } });
}

function tableTitle(t) {
  return para([new TextRun({ text: t, italics: true, font: FONT, size: BODY })],
              { spacing: { before: 0, after: 120, line: 480 } });
}

function note(parts) {
  return para(
    [new TextRun({ text: 'Note. ', italics: true, font: FONT, size: TBL })].concat(
      parts.map(s => new TextRun({
        text: s.t, italics: !!s.i, superScript: !!s.sup,
        font: FONT, size: TBL,
      }))),
    { spacing: { before: 120, after: 0, line: 240 } });
}

const blank = () => para([new TextRun({ text: '', font: FONT, size: BODY })]);

// ---------------------------------------------------------------- Table 1
// 8 columns. Usable width on Letter with 1 in margins = 12240 - 2880 = 9360 DXA.
const W1 = [1800, 700, 1000, 1200, 1100, 1100, 1250, 1210];
const sum1 = W1.reduce((a, b) => a + b, 0);

const t1 = new Table({
  columnWidths: W1,
  width: { size: sum1, type: WidthType.DXA },
  borders: {
    top: NONE, bottom: NONE, left: NONE, right: NONE,
    insideHorizontal: NONE, insideVertical: NONE,
  },
  rows: [
    // Header row 1: stubs span both header rows; two spanner heads.
    new TableRow({
      children: [
        cell('Modality', { width: W1[0], rowSpan: 2, align: AlignmentType.LEFT, top: true }),
        cell([{ t: 'n', i: true }], { width: W1[1], rowSpan: 2, top: true }),
        cell([{ t: 'R', i: true }, { t: '2', sup: true }], { width: W1[2], rowSpan: 2, top: true }),
        cell('MAE (years)', { width: W1[3], rowSpan: 2, top: true }),
        cell([{ t: 'r', i: true }, { t: ' (BAG, age)' }],
             { width: W1[4] + W1[5], span: 2, top: true, bottom: true }),
        cell([{ t: 'SD', i: true }, { t: ' of BAG (years)' }],
             { width: W1[6] + W1[7], span: 2, top: true, bottom: true }),
      ],
    }),
    // Header row 2: the column heads under each spanner.
    new TableRow({
      children: [
        cell('Before', { width: W1[4], bottom: true }),
        cell('After', { width: W1[5], bottom: true }),
        cell('Between persons', { width: W1[6], bottom: true }),
        cell('Within person', { width: W1[7], bottom: true }),
      ],
    }),
    new TableRow({
      children: [
        cell('Structural', { width: W1[0], align: AlignmentType.LEFT }),
        cell('1,380', { width: W1[1] }),
        cell('.82 (.03)', { width: W1[2] }),
        cell('5.35 (0.39)', { width: W1[3] }),
        cell('−.49', { width: W1[4] }),
        cell('.02', { width: W1[5] }),
        cell('7.42', { width: W1[6] }),
        cell('1.65', { width: W1[7] }),
      ],
    }),
    new TableRow({
      children: [
        cell('Functional', { width: W1[0], align: AlignmentType.LEFT, bottom: true }),
        cell('1,380', { width: W1[1], bottom: true }),
        cell('.69 (.04)', { width: W1[2], bottom: true }),
        cell('7.20 (0.39)', { width: W1[3], bottom: true }),
        cell('−.72', { width: W1[4], bottom: true }),
        cell('.08', { width: W1[5], bottom: true }),
        cell('9.94', { width: W1[6], bottom: true }),
        cell('3.40', { width: W1[7], bottom: true }),
      ],
    }),
  ],
});

const note1 = note([
  { t: 'N', i: true },
  { t: ' = 1,380. Accuracy was computed from the uncorrected out-of-sample predictions of the nested cross-validation and is reported as the mean over the 100 outer test folds (10 folds × 10 repeated partitions); values in parentheses are standard deviations across those folds. Variability across the 10 partitions was an order of magnitude smaller (' },
  { t: 'R', i: true }, { t: '2', sup: true },
  { t: ': ' }, { t: 'SD', i: true }, { t: ' = 0.003 and 0.004; MAE: ' }, { t: 'SD', i: true },
  { t: ' = 0.04 and 0.06 years, for the structural and functional models respectively). ' },
  { t: 'r', i: true },
  { t: ' (BAG, age) is the correlation between the brain-age gap and chronological age before and after the cross-fitted linear bias correction; values near zero after correction indicate that the age dependence of the gap has been removed. The between-persons ' },
  { t: 'SD', i: true },
  { t: ' is the dispersion of the corrected gap across participants and is the unit in which standardized coefficients in the subsequent analyses are expressed. The within-person ' },
  { t: 'SD', i: true },
  { t: ' is the mean standard deviation of a participant’s gap across the 10 partitions and indexes the stability of the measure. Structural = cortical and subcortical morphometry; Functional = tangent-space resting-state connectivity reduced by principal component analysis. The column heads Before and After refer to the bias correction. BAG = brain-age gap; MAE = mean absolute error.' },
]);

// ---------------------------------------------------------------- Table 2
const W2 = [2600, 1300, 1700, 1700, 1300];
const sum2 = W2.reduce((a, b) => a + b, 0);

const t2 = new Table({
  columnWidths: W2,
  width: { size: sum2, type: WidthType.DXA },
  borders: {
    top: NONE, bottom: NONE, left: NONE, right: NONE,
    insideHorizontal: NONE, insideVertical: NONE,
  },
  rows: [
    new TableRow({
      children: [
        cell('Modality', { width: W2[0], align: AlignmentType.LEFT, top: true, bottom: true }),
        cell([{ t: 'R', i: true }, { t: '2', sup: true }], { width: W2[1], top: true, bottom: true }),
        cell('MAE (years)', { width: W2[2], top: true, bottom: true }),
        cell('RMSE (years)', { width: W2[3], top: true, bottom: true }),
        cell([{ t: 'b', i: true }], { width: W2[4], top: true, bottom: true }),
      ],
    }),
    new TableRow({
      children: [
        cell('Structural', { width: W2[0], align: AlignmentType.LEFT }),
        cell('.78', { width: W2[1] }),
        cell('5.80', { width: W2[2] }),
        cell('7.42', { width: W2[3] }),
        cell('1.01', { width: W2[4] }),
      ],
    }),
    new TableRow({
      children: [
        cell('Functional', { width: W2[0], align: AlignmentType.LEFT, bottom: true }),
        cell('.61', { width: W2[1], bottom: true }),
        cell('7.90', { width: W2[2], bottom: true }),
        cell('9.94', { width: W2[3], bottom: true }),
        cell('1.05', { width: W2[4], bottom: true }),
      ],
    }),
  ],
});

const note2 = note([
  { t: 'N', i: true },
  { t: ' = 1,380. Values describe the accuracy of the bias-corrected predicted age and are pooled over all participants and all 10 partitions rather than averaged over folds, because the corrected predictions were not scored at the fold level. They are therefore not directly comparable with the fold-mean values in Table 1 and do not represent the predictive accuracy of the models, which Table 1 reports. ' },
  { t: 'b', i: true },
  { t: ' is the unstandardized slope of corrected predicted age regressed on chronological age; ' },
  { t: 'b', i: true },
  { t: ' ≈ 1 indicates that the regression-to-the-mean bias has been removed. The lower ' },
  { t: 'R', i: true }, { t: '2', sup: true },
  { t: ' relative to Table 1 is an expected consequence of the correction, which rescales the predictions to restore the slope against age at the cost of fit. Modalities are defined as in Table 1. MAE = mean absolute error; RMSE = root mean square error.' },
]);

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: BODY } } } },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840, orientation: PageOrientation.PORTRAIT },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: [
      tableNumber(1),
      tableTitle('Brain-Age Prediction Accuracy and Properties of the Brain-Age Gap'),
      t1,
      note1,
      blank(), blank(),
      tableNumber(2),
      tableTitle('Calibration Diagnostic: Accuracy of the Bias-Corrected Predicted Age'),
      t2,
      note2,
    ],
  }],
});

Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(process.argv[2], b);
  console.log('wrote ' + process.argv[2]);
});
