# n_epma_analysis

Repo containing data and analyses for nitrogen by EPMA project.

Notebooks:

- `01_buddingtonite_ANU.ipynb`: Fits wavescan, reads raw EPMA data, applies corrections, writes a CalcZAF input file, reads CalcZAF outputs to produce results
- `02_hyalophane_StA.ipynb`: Fits wavescan, reads raw EPMA data, applies corrections, writes a CalcZAF input file, reads CalcZAF outputs to produce results
- `03_basaltic_glasses_majors_wavescans.ipynb`: Compiles major element analyses for basaltic glasses, inspects and selects wavescans to fit, and performs fits for all glasses
- `04_basaltic_glasses_Nquant_D2872.ipynb`: Reads raw EPMA data, applies corrections, writes a CalcZAF input file, reads CalcZAF outputs to produce results


With reference to the paper, the figures and tables are listed below, pointing to the relevant notebooks to reproduce the data/figures.

Figures:

1. Figure 1: Buddingtonite + Tobelite count logs
2. Figure 2: Glass, Cs- and Pb-nitrate count logs
3. Figure 3: Hyalophane background scans
4. Figure 4: Four-point background correction method: `02_hyalophane_StA.ipynb`. This figure is based on a quantitative analysis of hyalophane so it is created in this notebook.
5. Figure 5: How APF was obtained for nitrides ???
6. Figure 6: How APF was obtained for glass and buddingtonite: `10_peak_shapes.ipynb`
7. Figure 7: How APF and uncertainty was calculated: `11_peak_shape_fits.ipynb`
8. Figure 8: APF results: `12_APF_relative_to_GaN.ipynb`
9. Figure 9: Peak positions
10. Figure 10: Measured vs reference values for the reference materials

Tables:

1. Table 1: List of reference materials
2. Table 2: APF results: `12_APF_relative_to_GaN.ipynb` and `12_APF_relative_to_BN.ipynb` 
3. Table 3: Buddingtonite quantitative results
4. Table 4: Hyalophane quantitative results
5. Table 5: Basaltic glasses - major/minor element quantitative results: `04_basaltic_glasses.ipynb`
6. Table 6: Basaltic glasses - nitrogen results
7. Table 7: Rhyolitic glasses - nitrogen results


## Peak shapes:

1. Run `10_peak_shapes.ipynb` to read and transform the data, produce figure 5
2. Run `11_peak_shape_fits.ipynb` to perform background and peak fits of original and MC simulated data, and produce figure 7
3. Run `12_APF_relative_to_GaN.ipynb` to calculate APF values, produce figure 8 and part of table 2.
4. Run `12_APF_relative_to_BN.ipynb` to calculate APF values relative to BN. This also produces an alternative version of figure 8, and additional numbers added to table 2.

APF values from table 2 are used in `data/_dictionaries/apf_values.csv` and `data/_dictionaries/apf_values_relative_to_BN.csv` for quantitative analysis.



## Calczaf usage

I was using CalcZAF v 13.2.4.

Note that some standards used in this analysis are not in the default Calczaf standard database. They were added manually. A custom standard database is included under `data/calczaf-standard-db-custom.mdb` (this is not added to the git repo however at present).

Be cautious about samples where N is specified as (NH4)2O: there's a bug in this version
of calczaf that means the input files are not correctly read in.

When you open the input file, it assigns 1 oxygen to each H cation, but there should be
zero oxygens for each H cation with this speciation.

Workaround: open the input file, manually edit the H valence selecting zero oxygens for each 1 cation. Close the file. This doesn't alter the file, but does alter the state
of CalcZAF for the remainder of your session so that it should read files correctly after this.