## DIFFERENT VERSIONS/FILES IN MY SUB-GRID MODEL

CO.txt - File used to read data for RT

debug_1.5.py - Running the model on the simulation; Lognormal-PDF and LTE
debug_1.5_plotting.ipynb - Results from debug_1.5.py; includes calculation of M_H2, L_CO, alpha_CO values for c_s_CO and delta_v values. Also includes plotting X_H2_bar, X_CO_bar, l_CO_bar/(n_H_mean*m_p) vs log(n_H_mean).

debug_1.6.py - Running the model on the simulation; Lognormal-PDF and non-LTE (just an approximation for non-LTE by using beta_nu to calculate T_exc)
debug_1.6_plotting.ipynb - Results from debug_1.6.py; includes calculation of M_H2, L_CO, alpha_CO values for c_s_CO and delta_v values. Also includes plotting X_H2_bar, X_CO_bar, l_CO_bar/(n_H_mean*m_p) vs log(n_H_mean).

debug_1.7.py - Running the model on the simulation; Power-law PDF and LTE
debug_1.7_plotting.ipynb - Results from debug_1.7.py; includes calculation of M_H2, L_CO, alpha_CO values for c_s_CO and delta_v values. Also includes plotting X_H2_bar, X_CO_bar, l_CO_bar/(n_H_mean*m_p) vs log(n_H_mean).

debug_1.8.py - Running the model on the simulation; Power-law PDF and non-LTE (just an approximation for non-LTE by using beta_nu to calculate T_exc)
debug_1.8_plotting.ipynb - Results from debug_1.8.py; includes calculation of M_H2, L_CO, alpha_CO values for c_s_CO and delta_v values. Also includes plotting X_H2_bar, X_CO_bar, l_CO_bar/(n_H_mean*m_p) vs log(n_H_mean).

The main files are debug_1.5 and debug_1.5_plotting. The plots and results (M_H2, L_CO, alpha_CO) in the report are from these files. The galaxy maps are only in debug_1.5_plotting, because galaxy maps for (Lognormal-PDF + LTE) only were needed. Although the same code can be used for different regimes across different plotting files listed above.

The results (M_H2, L_CO, alpha_CO) for the (Power-law PDF + LTE) in the report are from the file debug_1.7_plotting. Although plotting is done in the file, the plots were not included in the report. 


##Now the code behind these python scripts (you are probably looking for the above files, the following are just for looking at separate code for the sub-grid model):


sub_grid_model.ipynb - The code for the sub-grid model; Lognormal-PDF + LTE; Plotting also included.

pdf_evolution.ipynb - The code for computing the Power-law PDF; plotting also included.
