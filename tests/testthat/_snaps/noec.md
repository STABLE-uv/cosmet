# s11_criteria() uses the thresholds of Text S11

    Code
      print(s)
    Output
      <cosmet_criteria> Amador et al. (2024), Text S11: macroinvertebrates
      a. significant effect in a single direction
      b. minimum abundance: 3 individuals/sample
      c. monotonic means: yes (tolerance 0)
      d. no effect before the application: yes

# noec_table() lays out one row per taxon and one column per time

    Code
      print(noec)
    Output
      <meso_noec> from Williams' test (direction: both, alpha = 0.05)
      Criteria: none (significant effect only)
      Status: ok 2
      Effects detected: 1; selected: 1
      Selected NOECs (concentration unit of the data):
        taxon   -7  7
       Cloeon <NA> 1↓

