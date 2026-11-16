# check_design() counts replicates and reports problems

    Code
      print(chk)
    Output
      <meso_design_check>
      Replicates (mesocosms sampled) per treatment and time:
       treatment concentration 7 14
               0             0 3  3
               1             1 3  3
               3             3 3  3
              10            10 3  3
      
      Values: 24 recorded, 0 zero-filled, 0 missing, 0 not recorded
      No problems found.

# transform_response() computes ln(A x + 1) and records A

    Code
      print(d)
    Output
      <meso_data>
      12 mesocosms, 1 sampling time(s), 1 taxa or variables, 12 values
      Treatments: 0, 1, 3, 10; control: 0
      Times: 7
      Unit: not given
      Analysis scale: ln(2 x + 1)
      Values: 12 recorded, 0 zero-filled, 0 missing, 0 not recorded

