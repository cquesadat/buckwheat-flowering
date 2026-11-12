# Buckwheat flowering phenotyping and genetic analysis

Scripts and supporting data for image-based flowering phenotyping and genome-wide association analyses in common buckwheat (*Fagopyrum esculentum*).

## Related manuscript

This material accompanies the following preprint:

**bioRxiv preprint (2026):** [https://doi.org/10.64898/2026.11.03.763633](https://www.biorxiv.org/content/10.64898/2026.11.03.763633v1).

The preprint has not yet been peer reviewed. Once the peer-reviewed article is published, this README will be updated with the final publication citation and DOI.

## Download and setup

The analysis scripts and the accompanying data are distributed separately. Download both and **combine their contents in one project folder**, preserving the directory structure below. The enclosing folder name is arbitrary (`buckwheat-flowering` is only an example), and the scripts do not depend on it. Both downloads include the same top-level `README.md`, so it does not matter which copy is retained when merging.

```text
buckwheat-flowering/                  # Any folder name
├── README.md
├── A_Image_Segmentation_Modelling/
│   ├── README.md
│   ├── 0_Crop_Images_Experiment.sh
│   ├── ...                         # Analysis scripts and model code
│   └── Experiments/
│       ├── ETH23/
│       ├── ETH24/
│       ├── LfL23/
│       ├── LfL24/
│       └── ETH24L/
├── B_Flowering_dynamics_Modelling/
│   ├── README.md
│   ├── ...                         # Analysis scripts
│   ├── data/
│   ├── flowering_parameters/
│   └── BLUEs_BLUPs/
└── C_Genome_Wide_Analysis_Studies/
    ├── README.md
    ├── ...                         # Analysis scripts
    ├── data/
    └── outputs/
```

The A, B and C directories must be **siblings**. No path in the analysis code depends on the name of the enclosing folder or the names of the download sources. Run each pipeline from its own directory as shown in its README.

## Analyses

| Directory | Purpose | Main archived data/results |
| --- | --- | --- |
| `A_Image_Segmentation_Modelling/` | Mask R-CNN flower segmentation and image-derived flower cover | Original photographs, VIA polygon annotations and flower-cover time series |
| `B_Flowering_dynamics_Modelling/` | Flowering/abscission curve fitting, quality control and genotype estimates | Cleaned trait tables and genotype BLUEs/BLUPs |
| `C_Genome_Wide_Analysis_Studies/` | GWAS, gene annotations, marker-effect models and LD-supported comparisons | Genotype inputs, selected association, marker-effect and LD/concordance outputs |

### A. Image data

The five experiments are `ETH23`, `ETH24`, `LfL23`, `LfL24` and `ETH24L` (the latter is a separate second-crop trial). Each experiment has:

- `Design_<EXPERIMENT>.csv` — plot-to-accession mapping.
- `Original/Session_DD_MM_YYYY/` — original JPG images.
- `Predicted_images/Session_DD_MM_YYYY/` — reconstructed VIA-format flower polygons.
- `Results_Flowering_Perc_<EXPERIMENT>.csv` — retained measurements with columns `Plot`, `ID`, `Date`, `Total_Area`.

`Total_Area` is **flower-covered percentage of the centrally cropped image**, not the uncropped photograph. The provided polygon coordinates use the cropped-image reference frame. The archived percentage tables can be analysed without rerunning segmentation. Intermediate crops and individual quadrant predictions are not needed for those analyses.

At the time of packaging, the experiment-level filename audit matched **37,073 retained CSV records** to original JPGs and prediction JSONs. This was a file-coverage check; the segmentation polygons were validated separately during their reconstruction.

### B. Flowering dynamics and genotype estimates

The `data/` directory supplies sowing dates, manual exclusions and manually scored BBCH traits. `flowering_parameters/` holds fitted and cleaned flowering/abscission traits. `BLUEs_BLUPs/` includes the five experiments, `ETH_across_years`, `LfL_across_years`, `BBCH_2023`, `BBCH_2024` and `BBCH_across_years`.

Across-year BBCH outputs use `BLUPs_BLUEs_<trait>.csv` in the archived data. The supplied scripts use this naming convention and the GWAS runner also accepts `<trait>_BLUEs_BLUPs.csv` when the archived convention is absent.

The provided cleaned trait files may retain historical names; consult the B pipeline README for the supported stages and the software needed for refitting.

### C. Genetic analyses

The `data/` directory contains `GWAFFs_imputed.txt` (SNP-by-accession allele frequencies), `SNP_information.txt`, `annotated_h2.gff3` and `out.emapper.annotations.tsv`. The genomic marker matrix contains **allele frequencies**, rather than an already converted 0/1/2 genotype-dosage matrix.

The `outputs/` directory contains a **selected archive**, not every intermediate product. Key retained files include:

- `associations/all_Filter_GWAS_results_3PC_combined.csv`, `automated_unique_hits.csv`, `automated_snps_summary_stats.csv` and `snp_trait_gene_associations.csv`.
- `marker_effects/SNP_fit_summary_all_genes.csv`, `SNP_best_per_gene_trait_dataset.csv` and `SNP_best_overall_per_gene.csv`.
- `LD_decay/ld_decay_windows_by_chromosome.csv`.
- `concordance/snp_to_ld_block.csv`, `genome_ld_blocks.csv`, `ld_edges_used.csv`, `ld_windows_used.csv`, `direct_hits_SHARED_CATEGORY_RULE.csv`, `block_hit_pairs_DISTINCT_POSITIONS_ONLY.csv`, `block_hits_by_block_DISTINCT_POSITIONS_ONLY.csv` and `summary_DISTINCT_POSITION_LD_CORRECTED.csv`.

Intermediate files not listed above can be regenerated by the corresponding analysis scripts. **Do not assume every downstream script can start independently from only the archived outputs**: for example, the selected-SNP model script expects `outputs/associations/annotated_marker_trait_associations.csv`, which is generated by the preceding annotation step and is not part of the selected archive. Follow the scripts' input dependencies when starting at a later stage. The LD-block concordance analysis can use the supplied combined GWAS table, marker inputs and LD-decay windows.

## Running the pipelines

Follow the pipeline-specific READMEs for full commands and dependencies. In summary:

1. **A:** Use archived `Results_Flowering_Perc_*` tables, or run segmentation with the trained model and a separate working directory.
2. **B:** From `B_Flowering_dynamics_Modelling/`, run `0_Move_and_adjust_data.sh`, then the flowering-fitting, filtering and genotype-estimation scripts. Archived clean parameters and genotype estimates are already available for downstream analysis.
3. **C:** From `C_Genome_Wide_Analysis_Studies/`, run `1.Run_GWAS.sh` for new GWAS or start from `2_Aggregate_and_Filter_GWAS.R` if individual runs exist. With the archived combined association table, downstream scripts can be run according to their own prerequisites. The supplied LD-decay windows allow the LD-block concordance step without repeating LD-decay estimation.

Some combined genotype-estimate analyses require licensed **ASReml-R**; segmentation requires the specified legacy TensorFlow/TensorPack environment. See the corresponding pipeline READMEs.

## Interpretation

Automated flowering traits and manual BBCH traits describe related but non-identical phenotypes. Archived GWAS association summaries must be interpreted with the methods and significance criteria in the accompanying manuscript; an aggregated GAPIT `Filter` table is not by itself a list of BH-FDR-significant markers.

Please cite the bioRxiv preprint linked above when using these materials. This citation guidance will be updated when the peer-reviewed article is published.
