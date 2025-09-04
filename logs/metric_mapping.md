This analysis uses locked definitions for gene-level metrics:
- Sensitivity: TP_genes / total_TUSCO_genes.
- Precision: TP_genes / (TP_genes + FP_genes) at the locus level (FP_gene = locus with any non-reference detection and no TP).
- PDR: (TP_genes + FP_genes) / total_TUSCO_genes.
- FDR: FP_genes / (TP_genes + FP_genes).

Earlier project scripts computed some metrics at the transcript level (e.g., non-redundant precision) and PDR based on FSM/ISM only. We adopt the locus-level definitions above for all recomputed numbers.
