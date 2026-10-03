# *Escherichia coli* from diseased pigs: two lineage definitions compared

`input.csv` holds one row per isolate with readings, 481 in all: the ciprofloxacin call (1 when the MIC lies above 1 mg/L, the technical cut of `config_two_categories.yaml`), the sequence type, the phylogroup as the source publishes it and the phylogroup with ST23, ST88 and ST410 in phylogroup C (`phylogroup_clermont2013`). The source codes these three sequence types as B1; the quadruplex scheme of Clermont et al. (2013) places clonal complex 23 in phylogroup C. The comparison asks how much of the lineage share of the call each definition sees, and which part of the difference comes from relabelling the same isolates.

```bash
amr-clonalshare compare --input examples/ecoli_swine/input.csv --id-column isolate_id \
  --outcome ciprofloxacin_call --lineage-a sequence_type --lineage-b phylogroup \
  --output out/ecoli_comparison
amr-clonalshare compare --input examples/ecoli_swine/input.csv --id-column isolate_id \
  --outcome ciprofloxacin_call --lineage-a sequence_type --lineage-b phylogroup_clermont2013 \
  --output out/ecoli_comparison_clermont
```

The runs use the default seed (42) and budgets; their records are in `expected_comparison/` and `expected_comparison_clermont/`. `build_derived_tables.py` writes `input.csv` from the shipped tables, and `DATA_PROVENANCE.md` documents the collection.
