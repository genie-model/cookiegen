# **********************************************************************
# adapted from:	muffin.CB.p_worjh2.BASES.caoetal.SPIN
# base-config: 	16-level, Caoetal.2009 tracer selection
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# SPECIFY WHETHER RADIATIVE FORCING IS DERIVED FROM ATMOSPHERIC CO2	
# NOTE: 'y' requires that CO2 is being calculated in the first place ...
#       'n' can be used for a configuration either with or without explicitly calculated atmopsheric CO2
# e.g., ea_36=y -> raditive forcing is set by the calculated value of atmopsheric CO2
ea_36=y
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
# SPECIFY THE BIOLOGICAL EXPORT SCHEME
# e.g., bg_par_bio_prodopt='1N1T_PO4MM' -> as used in Cao et al. [2009]
bg_par_bio_prodopt='1N1T_PO4MM'
# SPECIFY BIOLOGICAL EXPORT PARAMETERS
# maximum rate of conversion of dissolved PO4 into organic matter by phytoplankton (mol kg-1 yr-1)
bg_par_bio_k0_PO4=8.9876e-006
# [PO4] M-M half-sat value (mol kg-1)
bg_par_bio_c0_PO4=8.9369e-007
#
# *** EXPORT RATIOS ****************************************************
#
# SPECIFY BULK CARBON EXPORT (AND/OR REDFIELD) RATIOS
# production fraction of dissolved organic matter
bg_par_bio_red_DOMfrac=0.66
# underlying export CaCO3 as a proportion of particulate organic matter (i.e., CaCO3/POC)
bg_par_bio_red_POC_CaCO3=0.0485
# exponent for modifier of CaCO3:POC export ratio
bg_par_bio_red_POC_CaCO3_pP=0.7440
#
# *** REMINERALIZATION *************************************************
#
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF DISSOLVED ORGANIC MATTER (DOM)
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=0.5
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF PARTICULATE ORGANIC MATTER (POM)
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=0.0557
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=589.9451
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF CaCO3
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.45
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1.8905e+003
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
#
# *** DATA SAVING ******************************************************
#
# (refer to Chapter: 'cookie model output')
# SPECIFY THE SCOPE OF DATA SAVING
# e.g., bg_par_data_save_level=10 -> includes full carbonate chemsitry output
bg_par_data_save_level=10
# SPECIFY THE INTERVAL AT WHICH TIME-SERIES DATA IS SAVED
# e.g., bg_par_data_save_sig_timeinterval=10.0 -> save time-series data every 10.0 years
# NOTE: set bg_par_data_save_sig_timeinterval=0.0 
#       in order to employ a file of save years to determine when time-series data is saved
bg_par_data_save_sig_timeinterval=1.0
# SPECIFY THE INTERVAL AT WHICH TIME-SLICE (SPATIAL) DATA IS SAVED
# e.g., bg_par_data_save_slice_timeinterval=1000.0 -> save time-slice data every 1000.0 years
# NOTE: set bg_par_data_save_slice_timeinterval=0.0 
#       in order to employ a file of save years to determine when time-slice data is saved
# NOTE: additionally set bg_par_infile_slice_name='save_timeslice_NONE.dat'
#       in order to force time-slice save at run end only
bg_par_data_save_slice_timeinterval=0.0
bg_par_infile_slice_name='save_timeslice_NONE.dat'
# SPECIFY USE OF MORE EXPLICIT (EXPERIMENT NAME BASED) netCDF FILENAMES
# e.g., bg_ctrl_ncout_expid_name=.false. -> assign default (short) netCDF filenames
bg_ctrl_ncout_expid_name=.false.
#
# *** FORCINGS *********************************************************
#
# SPEFICY ANY FORCINGS OF OCEAN/ATMOSPHERE/SEDIMENTS
# e.g., bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2" -> generic forcing of atmopsheric pCO2 and d13C
bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2"
# SCALE THE FORCING
bg_par_atm_force_scale_val_3=278.0E-06
bg_par_atm_force_scale_val_4=-6.5
#
# *** MISC **********************************************************
#
#
# *******************************************************************
# *** END ***********************************************************
# *******************************************************************

