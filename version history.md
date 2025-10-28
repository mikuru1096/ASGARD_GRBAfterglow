**v0.1**  Early-stage forward shock synchrotron radiation code.  
**v0.2**  Project initiation, named HE_pipeline, developed SSC module.  
**v0.3**  Developed $\gamma\gamma$ annihilation module.  
**v0.4**  Developed EBL absorption module.  
**v0.x**  Extensive code refactoring and bug fixes.  
**v1.0**  First version completed during the 2021 Spring Festival, applied to GRB 190829A.  
**v1.1**  Introduced Python interface.  
**v1.2**  Introduced $\chi^2$ module and extinction module.  
**v2.0**  Officially renamed ASGARD, stored in a private GitHub repository.  
**v2.1**  Improved electron spectrum calculation scheme, introduced CIP/WENO5/full-implicit algorithms.  
**v2.2**  Significant code optimization, introduced OpenMP.  
**v3.0**  August 2022: Abandoned the original interpolated electron spectrum approach, migrated to interpolated radiation spectrum scheme, improving program efficiency ($10^3$).  
**v3.1**  Introduced structured jet implementation based on Python multiprocessing.  
**v3.2**  Introduced forward and reverse shock implementation.  
**v3.2.1**  Fixed an issue where specific parameter combinations caused an excessively small initial afterglow time T00, leading to reverse shock dynamics calculation failure. (Feedback from Yun Wang)  
**v3.3**  Minor update (forward shock) based on a more realistic Compton Y parameter solution (Nakar).  
**v3.4**  Program efficiency optimization: merged synchrotron radiation module into electron spectrum calculation module.  

**v4.0**  Expanded IC cooling module to three methods (`dot_gam_e_ssc`, Nakar, Fan). The first two are nearly consistent, while Fan becomes unreliable in strong IC cooling regions.  
Speed ranking: Fan > Nakar > `dot_gam_e_ssc`; accuracy ranking is the opposite.  
**v4.0.1**  Fixed an issue in SED interpolation where flux values were incorrectly positioned at grid center i+1/2 instead of i-1/2, causing slightly elevated results (<1%). (Feedback from Jian-heZheng)  
**v4.0.2**  Added `dot_gam_e_ssa` module to calculate SSA-induced pile-up effects. (Beta version)
**v4.0.3**  Corrected the initial mass issue in Region 2 of reverse shock dynamics (Yan), which previously led to unphysical crossing timescales and incorrect light curve morphology. (Feedback from Yun Wang)  
**v4.0.4**  Resolved a historical issue with the Compton Y parameter, introduced in an uncertain version, which caused underestimation of `hat_gamma_e` in the Fan method. Cumulative updates.
