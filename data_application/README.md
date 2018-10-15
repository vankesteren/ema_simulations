### Data application: Epigenetics of childhood trauma and stress reactivity
Supplementary material to the manuscript _Exploratory Mediation Analysis with Many Potential Mediators_. For questions, contact E.-J. van Kesteren (e.vankesteren1@uu.nl).

#### How to run
It is possible to run each `.R` file separately as results are saved in the `results` folder.

0. Download the `*_sample_table.txt` files and the `samples.tsv` file from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-77445) into a folder called `data` in this folder.
1. Open `data_application.Rproj` in `RStudio`.
2. `preprocessing.R` loads the data and performs prefiltering.
3. `analysis.R` runs the CMF, HIMA and Filter algorithms.
4. `annotation.R` annotates the top sites from the CMF results.
5. `estimate.R` estimates the full mediation model for these sites.

#### Funding details
This work was supported by the Netherlands Organization for Scientific Research
(NWO) under Grant number `406.17.057`.
