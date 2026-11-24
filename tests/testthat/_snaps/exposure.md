# the exposure object prints a summary

    Code
      print(meso_exposure_small())
    Output
      <meso_exposure>
      6 mesocosms, 36 samples: 24 measured, 0 zero, 12 missing
      Treatments: 0 [control], 10, 100; unit ug/L
      Sampling times: 0, 7, 14 (before and after), 21, 35
      Applications at time 0, 14, each at the nominal concentration of the treatment

# the fit prints a summary

    Code
      print(fit_dissipation(e))
    Output
      <meso_dissipation>
      First-order dissipation: log-linear regression of ln C; one rate constant per application; one fit per treatment
      95% confidence intervals (Wald, parameter uncertainty only; those of k in $estimates); k per d, DT50 and DT90 in d
       treatment period n C0 (% nominal) k    DT50 DT50 CI      DT90 status
       10 ug/L   1      6 10 (100%)      0.1  6.93 [6.93, 6.93] 23   ok    
       10 ug/L   2      6 12 (120%)      0.05 13.9 [13.9, 13.9] 46.1 ok    
       100 ug/L  1      6 100 (100%)     0.1  6.93 [6.93, 6.93] 23   ok    
       100 ug/L  2      6 120 (120%)     0.05 13.9 [13.9, 13.9] 46.1 ok    
    Code
      print(fit_dissipation(e, method = "nonlinear", k_by = "series", by = "mesocosm"))
    Output
      <meso_dissipation>
      First-order dissipation: nonlinear least squares of C; one rate constant for all the applications; one fit per mesocosm
      95% confidence intervals (Wald, parameter uncertainty only; those of k in $estimates); k per d, DT50 and DT90 in d
       treatment mesocosm period n C0 (% nominal) k      DT50 DT50 CI    DT90 status
       10 ug/L   L1       1      3 9.02 (90%)     0.0646 10.7 [6.54, 30] 35.7 ok    
       10 ug/L   L1       2      3 12.5 (125%)    0.0646 10.7 [6.54, 30] 35.7 ok    
       10 ug/L   L2       1      3 9.02 (90%)     0.0646 10.7 [6.54, 30] 35.7 ok    
       10 ug/L   L2       2      3 12.5 (125%)    0.0646 10.7 [6.54, 30] 35.7 ok    
       100 ug/L  H1       1      3 90.2 (90%)     0.0646 10.7 [6.54, 30] 35.7 ok    
       100 ug/L  H1       2      3 125 (125%)     0.0646 10.7 [6.54, 30] 35.7 ok    
       100 ug/L  H2       1      3 90.2 (90%)     0.0646 10.7 [6.54, 30] 35.7 ok    
       100 ug/L  H2       2      3 125 (125%)     0.0646 10.7 [6.54, 30] 35.7 ok    

# each treatment follows its own schedule of applications

    Code
      print(e)
    Output
      <meso_exposure>
      5 mesocosms, 60 samples: 48 measured, 0 zero, 12 missing
      Treatments: control (0) [control], weekly (5), single (20); unit ug/L
      Sampling times: 0, 3, 7 (before and after), 10, 14 (before and after), 17, 21 (before and after), 24, 28
      Applications per treatment, time (nominal):
        weekly: 0 (5), 7 (5), 14 (5), 21 (5)
        single: 0 (20)

