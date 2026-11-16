# Step-down Williams test: reduced design against the variance of all groups

Written by `dev/validation/williams_variants.R` on 2026-11-16. Data `azoxy`, ln(Ax + 1) with the A of Amador et al. (2024), alpha = 0.05, both directions (p x 2), 99 999 draws, seed 123.

| Community | Habitat | Taxa x dates tested | Different NOECs | Detail (taxon, day: reduced -> all groups) |
|---|---|---|---|---|
| phytoplankton | eutrophic_lake | 118 | 9 | Achnanthes, 21: 20 -> <2; Amphora, 7: 20 -> <2; Cocconeis, 21: 20 -> <2; Coelastrum, 35: 20 -> 2; Coelastrum, 56: 20 -> 2; Cosmarium sp 2, 56: 20 -> <2; Cryptomonas sp, -7: <2 -> 2; Nitzschia spp, 56: 2 -> <2; Rhodomonas sp, 35: NA -> 2 |
| phytoplankton | water_spring | 135 | 10 | Achnanthes, 7: NA -> 20; Aphanocapsa, 35: 20 -> <2; Coccal < 5 um, -7: NA -> 20; Coelastrum, 35: 20 -> 2; Gomphonema, 56: 20 -> <2; Gyrosigma, -7: 20 -> <2; Oscillatoria, 56: 20 -> <2; Scenedesmus arcuatus, 21: 20 -> <2; Vorticella, 35: 20 -> <2; Vorticella, 7: NA -> 20 |
| zooplankton | eutrophic_lake | 124 | 2 | Pleuroxus_sp, 7: NA -> <2; Polyarthra, 35: NA -> 20 |
| zooplankton | water_spring | 106 | 15 | Brachionus angularis, 56: NA -> 2; Ceriodaphnia_sp, 7: 2 -> <2; Daphnia_sp, 21: 20 -> <2; Daphnia_sp, 56: 20 -> <2; Diaphanosoma, 56: NA -> 20; Euchlanis, 56: 20 -> <2; Filinia, 35: NA -> 20; Hydra, 21: 20 -> <2; Hydra, 35: 20 -> <2; Lecane luna, 56: NA -> 20; Lecane lunaris, 35: 2 -> <2; Lepadella, 7: NA -> 20; Macrothrix_sp, 56: NA -> 2; Trichocerca, 21: 20 -> <2; Trichocerca, 56: NA -> 20 |
| macroinvertebrates | eutrophic_lake | 77 | 2 | Colymbetes fuscus (adult), 35: 20 -> 2; Hydracarina, 7: NA -> 20 |
| macroinvertebrates | water_spring | 96 | 11 | Colymbetes fuscus larvae, 21: 20 -> <2; Colymbetes fuscus larvae, 35: NA -> <2; Dugastella valentina, 35: 2 -> <2; Dytiscidae (small larva), 56: NA -> 20; Eretes griseus (adult), 56: 20 -> <2; Lumbriculidae, 56: 20 -> <2; Micronecta, 21: 2 -> <2; Physella Acuta, 56: <2 -> 2; Theodoxus, 35: 20 -> <2; Zygoptera, 21: NA -> 20; Zygoptera, 56: 20 -> <2 |

`test_control(step_variance = "all")` against the independent implementation of this script (equal NOECs / taxa x dates tested):

- phytoplankton, eutrophic_lake: 118 of 118 NOECs equal
- phytoplankton, water_spring: 135 of 135 NOECs equal
- zooplankton, eutrophic_lake: 124 of 124 NOECs equal
- zooplankton, water_spring: 106 of 106 NOECs equal
- macroinvertebrates, eutrophic_lake: 77 of 77 NOECs equal
- macroinvertebrates, water_spring: 96 of 96 NOECs equal
