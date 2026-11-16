# print methods are stable

    Code
      print(tw)
    Output
      <meso_test> Williams' test, direction: both, alpha = 0.05
      Scale: ln(2 x + 1); 15 taxa x times
      Steps: variance and isotonic means of all groups; null distribution: 99,999 simulated draws, seed 123
      Status of the taxa x times: ok 14, absent 1
      Taxa x times with a significant difference from the control: 3
    Code
      print(derive_noec(tw, pre_times = -7, criteria = s11_criteria(
        "macroinvertebrates")))
    Output
      <meso_noec> from Williams' test (direction: both, alpha = 0.05)
      Criteria: Amador et al. (2024), Text S11: macroinvertebrates
      Status: ok 14, absent 1
      Effects detected: 3; selected: 2
      Selected NOECs (concentration unit of the data):
             taxon   -7    7   21   35   56
       Chironomini <NA> <NA> <NA> <NA>  20↓
        Micronecta <NA> <NA>  20↓ <NA> <NA>

