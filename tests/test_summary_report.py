#!/usr/bin/env python3
"""Tests for the standalone Step 27 concise summary report."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO_DIR = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_DIR / "bin" / "Step_27_generate_summary_report.py"
SPEC = importlib.util.spec_from_file_location("evscope_summary_report", SCRIPT_PATH)
SUMMARY = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(SUMMARY)


class SummaryReportTests(unittest.TestCase):
    def test_paired_total_rna_summary_contains_requested_metrics(self):
        metrics = {
            "Total Raw Read Pairs (Fragments)": "42650000",
            "Total Trimmed Reads (R1+R2)": "78361898",
            "Number of Reads after UMI-deduplication (STAR Input)": "1030000",
            "Number of Uniquely Mapped Reads (STAR)": "440000",
            "Number of Multi-mapped Reads (STAR)": "485000",
            "Number of Reads Mapped to Escherichia coli (BBSplit)": "9081000",
            "Percentage of Trimmed Reads Mapped to E. coli (BBSplit)": "11.59",
            "Number of Reads Mapped to Mycoplasma (BBSplit)": "7766000",
            "Percentage of Trimmed Reads Mapped to Mycoplasma (BBSplit)": "9.91",
            "Number of Reads Mapped to Intronic Regions": "510000",
            "Percentage of Reads Mapped to Intronic Regions": "55.22",
            "Number of Reads Mapped to Exonic Regions": "58000",
            "Percentage of Reads Mapped to Exonic Regions": "6.30",
            "Number of Reads Mapped to Intergenic Regions": "278000",
            "Percentage of Reads Mapped to Intergenic Regions": "30.05",
            "Percentage of Mapped Reads on Forward Strand": "49.90",
            "Percentage of Mapped Reads on Reverse Strand": "50.10",
        }
        genes = [(f"GENE{i}", float(1100 - i)) for i in range(1, 11)]
        result = SUMMARY.build_summary(
            "sample_A", metrics, genes, "smarter_total_rna", "unstrand", True
        )

        self.assertIn("42.65 million read pairs", result)
        self.assertIn("39.18 million read pairs remaining after trimming", result)
        self.assertIn("1.03 million complete read pairs", result)
        self.assertIn("925,000 (89.81%) mapped to hg38", result)
        self.assertIn("47.57% were unique and 52.43% were multi-mapped", result)
        self.assertIn("11.59% of trimmed read pairs", result)
        self.assertIn("9.91% of trimmed read pairs", result)
        self.assertIn("55.22%", result)
        self.assertIn("6.30%", result)
        self.assertIn("30.05%", result)
        self.assertIn("50.10% of strand-informative alignments", result)
        self.assertIn("top 10 genes by read count", result)
        self.assertIn("GENE10", result)

    def test_expression_parser_excludes_circrna_and_aggregates_symbols(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            matrix = Path(tmpdir) / "expression.tsv"
            matrix.write_text(
                "GeneID\tGeneSymbol\tGeneType\tReadCounts\tNorm_Expr\n"
                "ENSG1\tTTN\tprotein_coding\t100\t1\n"
                "ENSG2\tTTN\tprotein_coding\t25\t1\n"
                "circ1\tcirc-X\tcircRNAs\t1000\t1\n"
                "ENSG3\tRERE\tprotein_coding\t80\t1\n",
                encoding="utf-8",
            )
            genes = SUMMARY.top_genes(str(matrix))

        self.assertEqual(genes, [("TTN", 125.0), ("RERE", 80.0)])

    def test_smrna_summary_reports_no_dedup_and_forward_orientation(self):
        metrics = {
            "Total Raw Reads (R1+R2)": "1000",
            "Total Trimmed Reads (R1+R2)": "900",
            "Number of Reads in Refined STAR Input (No UMI Deduplication)": "900",
            "Number of Uniquely Mapped Reads (STAR)": "700",
            "Number of Multi-mapped Reads (STAR)": "100",
            "Percentage of Mapped Reads on Forward Strand": "72",
            "Percentage of Mapped Reads on Reverse Strand": "18",
        }
        result = SUMMARY.build_summary(
            "smrna_sample", metrics, [], "smarter_smrna", "forward", False
        )

        self.assertIn("No UMI deduplication was applied", result)
        self.assertIn("80.00% of strand-informative alignments", result)
        self.assertIn("expected forward/sense orientation", result)

    def test_html_contains_one_summary_paragraph(self):
        report = SUMMARY.render_html("sample<&>", "A short & safe summary.")
        self.assertEqual(report.count('<p class="summary">'), 1)
        self.assertIn("sample&lt;&amp;&gt;", report)
        self.assertIn("A short &amp; safe summary.", report)


if __name__ == "__main__":
    unittest.main()
