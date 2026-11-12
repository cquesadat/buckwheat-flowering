# B - Flowering dynamics and genotype estimates

This pipeline fits flowering and abscission curves to image-derived flower cover, filters plot-level estimates, and calculates genotype BLUEs/BLUPs. It also processes manually scored BBCH traits.

## Inputs

- `../A_Image_Segmentation_Modelling/Experiments/<EXPERIMENT>/Results_Flowering_Perc_<EXPERIMENT>.csv` — `Plot,ID,Date,Total_Area`
- `../A_Image_Segmentation_Modelling/Experiments/<EXPERIMENT>/Design_<EXPERIMENT>.csv` — field design
- `data/sowing_dates.csv` — experiment and sowing date (`Trial,Sowing_Date`)
- `data/manual_exclusions.csv` — exclusions (`experiment,scope,field,value`)
- `data/parameter_bbch_2023.csv` and `data/parameter_bbch_2024.csv` — manual phenology traits

`<EXPERIMENT>` is one of `ETH23`, `ETH24`, `LfL23`, `LfL24` or `ETH24L`.

## Run

From `B_Flowering_dynamics_Modelling/`:

```bash
bash 0_Move_and_adjust_data.sh
Rscript 1_Flowering_Model_Fitting.R
Rscript 1.5_Filter_and_Convert_to_DAS.R
Rscript 2_BLUEs_BLUPs_all_experiments_and_BBCH.R
```

The first step stages copies of A's measurements and designs in `data/` and computes days after sowing (DAS). The next two fit and quality-filter plot-level flowering and abscission parameters. The final script estimates genotype BLUEs/BLUPs by experiment and across years, including the BBCH traits.

## Outputs

- `flowering_parameters/` — fitted parameters, fit metrics, cleaned DAS traits and quality-control tables
- `BLUEs_BLUPs/ETH23/`, `ETH24/`, `LfL23/`, `LfL24/`, `ETH24L/` — experiment-level estimates
- `BLUEs_BLUPs/ETH_across_years/`, `LfL_across_years/` — combined estimates
- `BLUEs_BLUPs/BBCH_2023/`, `BBCH_2024/`, `BBCH_across_years/` — manual phenology estimates

Across-year BBCH files are written as `BLUPs_BLUEs_<trait>.csv` to match the archived tables and Pipeline C. Other groups use `<trait>_BLUEs_BLUPs.csv`.

The supplied `flowering_parameters/` and `BLUEs_BLUPs/` directories contain archived results and can be used without refitting the models. The five archived cleaned parameter files retain names such as `parameters_clean_das_2023_ETH.csv`; the spatial-analysis script recognizes these alongside newly generated `parameters_<EXPERIMENT>_clean_das.csv` files.

Requires R and packages used by the scripts, including `dplyr`, `lubridate`, `SpATS` and licensed `ASReml-R` for combined analyses.