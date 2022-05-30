# Gna-EW-large-ensemble
Model code and analysis/figuremaking scripts for Gnanadesikan+WASP EW work
********************************************************************************
********************************************************************************
The model code has been run with GNA_WASP_F2_emissions with the following physical parameters altered:
1. Background emissions - CO2 emissions based off of Meinshausen et al 2011 for RCP  2.6, 4.5, 6.0, or 8.5. Yearly emissions from 1860-2100 are read in from XXX.XXXX.
2. Enhanced weathering - either 0 PgC/yr drawdown (control), 2 PgC/yr drawdown for years 2000-2100 (EW-FAST), or 1 PgC/yr drawdown for years 2000-2200 (EW-SLOW). EW-SLOW is determined by putting in the same drawdown rate for EW-FAST but setting the half_time value to 1.
3. Climate feedback parameter (lambda)
4. Southern Ocean wind stress (tau)
5. Southern Ocean subduction fraction (delta)

When terrestrial carbon dynamics are turned ON, the values for the following are set:
a. The CO2 fertilization coefficient (gamma_co2) to 0.5
b. The linear dependence of NPP with temperature (dNPP_dT) to -2 
c. The linear dependence of soil respiration with temperature (dtau_dT) to -0.5

When terrestrial carbon dynamics are turned OFF, the values for a-c are set to 0.
********************************************************************************
********************************************************************************
I save the model output in the following filename style: numbers following d, t, and l are values for delta, tau, and lambda (times 10 or 100 or 1000 as I did not name everything the same, apologies!)

Terrestrial dynamics OFF:
g2_rcp85_ctrl_d30_t120_l06.mat --> control run
g2_rcp85_2w2000_d30_t120_l06.mat --> EW-FAST run
g2_rcp85_2w2000_d30_t120_l06_half_rate.mat --> EW-SLOW run

Terrestrial dynamics ON:
terr_r85l08t08d05_ctrl.mat --> control run
terr_r85l08t08d05_ewf.mat --> EW-FAST run
terr_r85l08t08d05_ews.mat --> EW-SLOW run
