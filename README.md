# n_epma_analysis

Repo containing data and analyses for nitrogen by EPMA project.

With reference to the paper, the figures and tables are listed below, pointing to the relevant notebooks to reproduce the data/figures.

Figures:

1. Figure 1: Buddingtonite + Tobelite count logs
2. Figure 2: Glass, Cs- and Pb-nitrate count logs
3. Figure 3: Hyalophane background scans
4. Figure 4: Four-point background correction method
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
5. Table 5: Basaltic glasses - major/minor element quantitative results
6. Table 6: Basaltic glasses - nitrogen results
7. Table 7: Rhyolitic glasses - nitrogen results


## Peak shapes:

1. Run `10_peak_shapes.ipynb` to read and transform the data, produce figure 5
2. Run `11_peak_shape_fits.ipynb` to perform background and peak fits of original and MC simulated data, and produce figure 7
3. Run `12_APF_relative_to_GaN.ipynb` to calculate APF values, produce figure 8 and part of table 2.
4. Run `12_APF_relative_to_BN.ipynb` to calculate APF values relative to BN. This also produces an alternative version of figure 8, and additional numbers added to table 2.

APF values from table 2 are used in `data/_dictionaries/apf_values.csv` and `data/_dictionaries/apf_values_relative_to_BN.csv` for quantitative analysis.


## 


