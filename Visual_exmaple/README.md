# Visual example

This folder contains an anonymized EVscope report from a paired-end SMARTer smRNA analysis. It illustrates the report layout, embedded FastQC and Krona reports, interactive QC tables, RNA composition figures, and coverage plots.

| File | Description |
| --- | --- |
| [EVscope_example_report.html](EVscope_example_report.html) | Self-contained interactive HTML report, approximately 28 MB. |
| [report_preview.png](report_preview.png) | Screenshot of the strandness and splice-junction QC section. |
| [validation.json](validation.json) | Checks covering identifiers, retained table values, figures, embedded reports, and offline browser behavior. |

To view the complete report, open the HTML file on GitHub, click **Download raw file**, and open the downloaded file in a browser. A repository clone also contains everything needed to view it offline.

Names, project/sample identifiers, sequencing file identifiers, private filesystem paths, and run dates/times have been removed or replaced, including identifiers inside embedded reports, figure labels, and downloadable images. The example uses the sample identifier `Example_001`.

Original analysis values, gene and species names, sequence content, public software versions, and license information are retained. This is an example with identifying labels removed; it contains real analysis results. The original report was left unchanged, and no original-to-example identifier mapping is included.

![Strandness and splice-junction QC in the anonymized EVscope report](report_preview.png)
