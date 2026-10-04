# Same-host comparison, part E: field edge at 100 and 150 MeV (df2fabbd2566fc38d36dfd32daecc3147fa2403c)

Same host, same thread count, and the same fix in both codes. Margins are those of the 200 MeV field (#31). One field (15 x 15 cm) in a uniform medium, no range shifter. Every endpoint was computed by the analysis from hash-verified dose files.

## A-port against A-up (confirmatory)

- Joint claim (all 34 endpoints equivalent): **ESTABLISHED**
- Secondary joint claim (every endpoint other than rows where every run of both arms is zero is equivalent after Holm over all 34): **ESTABLISHED**; rows excluded: none
- Equivalent: 34 of 34 unadjusted, 34 after Holm; not established: 0

| endpoint | arm mean | reference mean | estimate | 90% interval | margin | outcome | Holm p |
|---|---|---|---|---|---|---|---|
| E100/lateral_39_5 | 4.2525 | 4.2727 | -0.0202 | [-0.0365, -0.00389] | [-0.5000, 0.5000] (difference) | equivalent | 1.34e-15 |
| E100/lateral_39_10 | 2.1167 | 2.0992 | 0.0175 | [0.00766, 0.0274] | [-0.2000, 0.2000] (difference) | equivalent | 7.87e-11 |
| E100/lateral_39_20 | 1.1045 | 1.1033 | 0.00126 | [-0.00915, 0.0117] | [-0.2000, 0.2000] (difference) | equivalent | 8.83e-14 |
| E100/lateral_39_30 | 0.6784 | 0.6774 | 0.000963 | [-0.0115, 0.0135] | [-0.2000, 0.2000] (difference) | equivalent | 1e-12 |
| E100/lateral_39_50 | 0.2472 | 0.2494 | -0.00227 | [-0.00553, 0.000999] | [-0.2000, 0.2000] (difference) | equivalent | 1.17e-17 |
| E100/lateral_39_70 | 0.0831 | 0.0829 | 0.000183 | [-0.00253, 0.00289] | [-0.2000, 0.2000] (difference) | equivalent | 6.46e-19 |
| E100/lateral_61_5 | 5.9305 | 5.9578 | -0.0274 | [-0.0444, -0.0103] | [-0.5000, 0.5000] (difference) | equivalent | 5.92e-14 |
| E100/lateral_61_10 | 2.6533 | 2.6512 | 0.00213 | [-0.0120, 0.0163] | [-0.2000, 0.2000] (difference) | equivalent | 6.65e-11 |
| E100/lateral_61_20 | 1.1911 | 1.1950 | -0.0039 | [-0.0156, 0.00778] | [-0.2000, 0.2000] (difference) | equivalent | 4.13e-13 |
| E100/lateral_61_30 | 0.7120 | 0.7122 | -0.000257 | [-0.00817, 0.00765] | [-0.2000, 0.2000] (difference) | equivalent | 1.82e-15 |
| E100/lateral_61_50 | 0.2645 | 0.2668 | -0.00237 | [-0.00715, 0.0024] | [-0.2000, 0.2000] (difference) | equivalent | 2.57e-18 |
| E100/lateral_61_70 | 0.0889 | 0.0899 | -0.001 | [-0.0043, 0.0023] | [-0.2000, 0.2000] (difference) | equivalent | 5.58e-17 |
| E100/cax_39 | 50940.0132 | 50990.0732 | 0.9990 | [0.9982, 0.9998] | [0.9950, 1.0050] (ratio) | equivalent | 1.74e-06 |
| E100/cax_61 | 70248.2568 | 70246.0508 | 1.0000 | [0.9989, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 3.94e-06 |
| E100/cax_max | 173554.7012 | 173560.8730 | 1.0000 | [0.9990, 1.0009] | [0.9950, 1.0050] (ratio) | equivalent | 7.7e-07 |
| E100/r80_mm | 75.7974 | 75.7970 | 0.000351 | [-0.000734, 0.00144] | [-0.5000, 0.5000] (difference) | equivalent | 2.01e-32 |
| E100/r20_mm | 78.3156 | 78.3156 | -1.87e-05 | [-0.00119, 0.00116] | [-0.5000, 0.5000] (difference) | equivalent | 2.5e-32 |
| E150/lateral_79_5 | 4.3813 | 4.3628 | 0.0184 | [-0.00474, 0.0416] | [-0.5000, 0.5000] (difference) | equivalent | 8.15e-14 |
| E150/lateral_79_10 | 2.6138 | 2.6132 | 0.000654 | [-0.0187, 0.0200] | [-0.2000, 0.2000] (difference) | equivalent | 1.9e-10 |
| E150/lateral_79_20 | 1.0437 | 1.0425 | 0.00121 | [-0.0098, 0.0122] | [-0.2000, 0.2000] (difference) | equivalent | 8.74e-13 |
| E150/lateral_79_30 | 0.4115 | 0.4088 | 0.00267 | [-0.00271, 0.00804] | [-0.2000, 0.2000] (difference) | equivalent | 1.5e-17 |
| E150/lateral_79_50 | 0.0541 | 0.0549 | -0.000789 | [-0.00273, 0.00115] | [-0.2000, 0.2000] (difference) | equivalent | 9.33e-22 |
| E150/lateral_79_70 | 0.00668 | 0.00661 | 6.98e-05 | [-0.00119, 0.00133] | [-0.2000, 0.2000] (difference) | equivalent | 1.02e-22 |
| E150/lateral_125_5 | 8.4315 | 8.4372 | -0.00579 | [-0.0358, 0.0242] | [-0.5000, 0.5000] (difference) | equivalent | 3.03e-11 |
| E150/lateral_125_10 | 3.6008 | 3.6213 | -0.0206 | [-0.0422, 0.000997] | [-0.2000, 0.2000] (difference) | equivalent | 8.66e-08 |
| E150/lateral_125_20 | 1.2799 | 1.2820 | -0.0021 | [-0.0137, 0.00952] | [-0.2000, 0.2000] (difference) | equivalent | 8.74e-13 |
| E150/lateral_125_30 | 0.4764 | 0.4801 | -0.00378 | [-0.0103, 0.00269] | [-0.2000, 0.2000] (difference) | equivalent | 2.97e-13 |
| E150/lateral_125_50 | 0.0573 | 0.0580 | -0.000756 | [-0.00218, 0.000671] | [-0.2000, 0.2000] (difference) | equivalent | 5.24e-23 |
| E150/lateral_125_70 | 0.0045 | 0.00462 | -0.000121 | [-0.000719, 0.000476] | [-0.2000, 0.2000] (difference) | equivalent | 3.65e-30 |
| E150/cax_79 | 36484.0024 | 36484.1099 | 1.0000 | [0.9989, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 3.94e-06 |
| E150/cax_125 | 46456.1875 | 46432.0713 | 1.0005 | [0.9995, 1.0015] | [0.9950, 1.0050] (ratio) | equivalent | 3.94e-06 |
| E150/cax_max | 112044.8701 | 111911.9102 | 1.0012 | [1.0000, 1.0024] | [0.9950, 1.0050] (ratio) | equivalent | 7.74e-05 |
| E150/r80_mm | 156.1147 | 156.1177 | -0.00293 | [-0.00603, 0.000173] | [-0.5000, 0.5000] (difference) | equivalent | 3.14e-25 |
| E150/r20_mm | 159.6966 | 159.7020 | -0.00538 | [-0.0107, -0.000107] | [-0.5000, 0.5000] (difference) | equivalent | 6.52e-21 |

Descriptive only. Ratio of arm means (arm / reference) at 50 and 70 mm outside the edge, 95% interval:

| endpoint | arm mean | reference mean | ratio | 95% interval | note |
|---|---|---|---|---|---|
| E100/lateral_39_50 | 0.2472 | 0.2494 | 0.9909 | [0.9751, 1.0069] |  |
| E100/lateral_39_70 | 0.0831 | 0.0829 | 1.0022 | [0.9630, 1.0430] |  |
| E100/lateral_61_50 | 0.2645 | 0.2668 | 0.9911 | [0.9697, 1.0130] |  |
| E100/lateral_61_70 | 0.0889 | 0.0899 | 0.9889 | [0.9452, 1.0345] |  |
| E150/lateral_79_50 | 0.0541 | 0.0549 | 0.9856 | [0.9435, 1.0296] |  |
| E150/lateral_79_70 | 0.00668 | 0.00661 | 1.0106 | [0.8016, 1.2741] |  |
| E150/lateral_125_50 | 0.0573 | 0.0580 | 0.9870 | [0.9575, 1.0174] |  |
| E150/lateral_125_70 | 0.0045 | 0.00462 | 0.9737 | [0.8298, 1.1427] |  |

Descriptive only. Beside each maximum: the R80 difference and the ratio of the three-row mean around each run's maximum row. No cause is attributed.

| energy | maximum: outcome | maximum: ratio, 90% | R80 difference (mm), 90% | three-row mean: ratio, 95% |
|---|---|---|---|---|
| 100 MeV | equivalent | 0.99996 [0.99902, 1.00091] | 0.000351 [-0.000734, 0.00144] | 0.99972 [0.99875, 1.00070]  |
| 150 MeV | equivalent | 1.00119 [0.99999, 1.00239] | -0.00293 [-0.00603, 0.000173] | 1.00086 [0.99942, 1.00230]  |

## B-pgcc against B-up (confirmatory)

- Joint claim (all 34 endpoints equivalent): **ESTABLISHED**
- Secondary joint claim (every endpoint other than rows where every run of both arms is zero is equivalent after Holm over all 34): **ESTABLISHED**; rows excluded: none
- Equivalent: 34 of 34 unadjusted, 34 after Holm; not established: 0

| endpoint | arm mean | reference mean | estimate | 90% interval | margin | outcome | Holm p |
|---|---|---|---|---|---|---|---|
| E100/lateral_39_5 | 4.2623 | 4.2612 | 0.00114 | [-0.0194, 0.0217] | [-0.5000, 0.5000] (difference) | equivalent | 6.02e-15 |
| E100/lateral_39_10 | 2.1178 | 2.1133 | 0.0045 | [-0.00948, 0.0185] | [-0.2000, 0.2000] (difference) | equivalent | 4.08e-12 |
| E100/lateral_39_20 | 1.1000 | 1.0994 | 0.000554 | [-0.00981, 0.0109] | [-0.2000, 0.2000] (difference) | equivalent | 8.22e-13 |
| E100/lateral_39_30 | 0.6814 | 0.6847 | -0.00334 | [-0.00857, 0.00189] | [-0.2000, 0.2000] (difference) | equivalent | 1.88e-17 |
| E100/lateral_39_50 | 0.2471 | 0.2479 | -0.000787 | [-0.00466, 0.00308] | [-0.2000, 0.2000] (difference) | equivalent | 6.1e-14 |
| E100/lateral_39_70 | 0.0852 | 0.0827 | 0.00255 | [-0.000734, 0.00583] | [-0.2000, 0.2000] (difference) | equivalent | 3.66e-20 |
| E100/lateral_61_5 | 5.9437 | 5.9403 | 0.00339 | [-0.0226, 0.0294] | [-0.5000, 0.5000] (difference) | equivalent | 1.75e-13 |
| E100/lateral_61_10 | 2.6589 | 2.6531 | 0.00589 | [-0.0093, 0.0211] | [-0.2000, 0.2000] (difference) | equivalent | 1.09e-11 |
| E100/lateral_61_20 | 1.1869 | 1.1944 | -0.00747 | [-0.0164, 0.00144] | [-0.2000, 0.2000] (difference) | equivalent | 2.47e-14 |
| E100/lateral_61_30 | 0.7126 | 0.7104 | 0.00224 | [-0.00518, 0.00966] | [-0.2000, 0.2000] (difference) | equivalent | 2.18e-14 |
| E100/lateral_61_50 | 0.2663 | 0.2618 | 0.0045 | [0.000988, 0.00801] | [-0.2000, 0.2000] (difference) | equivalent | 1.43e-18 |
| E100/lateral_61_70 | 0.0893 | 0.0879 | 0.0014 | [-0.00144, 0.00424] | [-0.2000, 0.2000] (difference) | equivalent | 2.19e-21 |
| E100/cax_39 | 50966.1626 | 50958.6147 | 1.0001 | [0.9991, 1.0012] | [0.9950, 1.0050] (ratio) | equivalent | 3.16e-06 |
| E100/cax_61 | 70255.7920 | 70252.2012 | 1.0001 | [0.9991, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 5.06e-06 |
| E100/cax_max | 173556.9570 | 173549.5547 | 1.0000 | [0.9990, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 4.75e-06 |
| E100/r80_mm | 75.7975 | 75.7976 | -0.000106 | [-0.00112, 0.000912] | [-0.5000, 0.5000] (difference) | equivalent | 7.32e-33 |
| E100/r20_mm | 78.3160 | 78.3159 | 8.04e-05 | [-0.00107, 0.00123] | [-0.5000, 0.5000] (difference) | equivalent | 1.5e-32 |
| E150/lateral_79_5 | 4.3548 | 4.3702 | -0.0154 | [-0.0368, 0.00604] | [-0.5000, 0.5000] (difference) | equivalent | 3.99e-13 |
| E150/lateral_79_10 | 2.6187 | 2.6228 | -0.00402 | [-0.0226, 0.0146] | [-0.2000, 0.2000] (difference) | equivalent | 1.2e-10 |
| E150/lateral_79_20 | 1.0475 | 1.0471 | 0.00039 | [-0.0118, 0.0126] | [-0.2000, 0.2000] (difference) | equivalent | 2.64e-11 |
| E150/lateral_79_30 | 0.4046 | 0.4059 | -0.00131 | [-0.00847, 0.00585] | [-0.2000, 0.2000] (difference) | equivalent | 1.51e-14 |
| E150/lateral_79_50 | 0.0538 | 0.0550 | -0.00118 | [-0.00409, 0.00174] | [-0.2000, 0.2000] (difference) | equivalent | 7.91e-19 |
| E150/lateral_79_70 | 0.00636 | 0.00624 | 0.000121 | [-0.000673, 0.000915] | [-0.2000, 0.2000] (difference) | equivalent | 1.09e-24 |
| E150/lateral_125_5 | 8.4438 | 8.4325 | 0.0113 | [-0.0141, 0.0367] | [-0.5000, 0.5000] (difference) | equivalent | 8.64e-12 |
| E150/lateral_125_10 | 3.6145 | 3.6261 | -0.0116 | [-0.0319, 0.00879] | [-0.2000, 0.2000] (difference) | equivalent | 8.18e-10 |
| E150/lateral_125_20 | 1.2870 | 1.2820 | 0.00505 | [-0.00507, 0.0152] | [-0.2000, 0.2000] (difference) | equivalent | 4.34e-13 |
| E150/lateral_125_30 | 0.4767 | 0.4791 | -0.00245 | [-0.0105, 0.00563] | [-0.2000, 0.2000] (difference) | equivalent | 1.97e-12 |
| E150/lateral_125_50 | 0.0572 | 0.0574 | -0.000222 | [-0.00195, 0.0015] | [-0.2000, 0.2000] (difference) | equivalent | 1.44e-23 |
| E150/lateral_125_70 | 0.00457 | 0.00452 | 5.02e-05 | [-0.000681, 0.000782] | [-0.2000, 0.2000] (difference) | equivalent | 3.13e-22 |
| E150/cax_79 | 36444.1035 | 36458.0591 | 0.9996 | [0.9985, 1.0008] | [0.9950, 1.0050] (ratio) | equivalent | 1.77e-05 |
| E150/cax_125 | 46440.3403 | 46437.5942 | 1.0001 | [0.9986, 1.0015] | [0.9950, 1.0050] (ratio) | equivalent | 2.36e-05 |
| E150/cax_max | 111967.8330 | 111967.8193 | 1.0000 | [0.9989, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 4.75e-06 |
| E150/r80_mm | 156.1136 | 156.1132 | 0.000425 | [-0.00264, 0.00349] | [-0.5000, 0.5000] (difference) | equivalent | 4.71e-25 |
| E150/r20_mm | 159.6991 | 159.6978 | 0.00129 | [-0.00246, 0.00504] | [-0.5000, 0.5000] (difference) | equivalent | 2.96e-19 |

Descriptive only. Ratio of arm means (arm / reference) at 50 and 70 mm outside the edge, 95% interval:

| endpoint | arm mean | reference mean | ratio | 95% interval | note |
|---|---|---|---|---|---|
| E100/lateral_39_50 | 0.2471 | 0.2479 | 0.9968 | [0.9778, 1.0163] |  |
| E100/lateral_39_70 | 0.0852 | 0.0827 | 1.0308 | [0.9830, 1.0810] |  |
| E100/lateral_61_50 | 0.2663 | 0.2618 | 1.0172 | [1.0007, 1.0339] |  |
| E100/lateral_61_70 | 0.0893 | 0.0879 | 1.0159 | [0.9770, 1.0564] |  |
| E150/lateral_79_50 | 0.0538 | 0.0550 | 0.9786 | [0.9162, 1.0452] |  |
| E150/lateral_79_70 | 0.00636 | 0.00624 | 1.0194 | [0.8729, 1.1906] |  |
| E150/lateral_125_50 | 0.0572 | 0.0574 | 0.9961 | [0.9602, 1.0334] |  |
| E150/lateral_125_70 | 0.00457 | 0.00452 | 1.0111 | [0.8304, 1.2312] |  |

Descriptive only. Beside each maximum: the R80 difference and the ratio of the three-row mean around each run's maximum row. No cause is attributed.

| energy | maximum: outcome | maximum: ratio, 90% | R80 difference (mm), 90% | three-row mean: ratio, 95% |
|---|---|---|---|---|
| 100 MeV | equivalent | 1.00004 [0.99895, 1.00114] | -0.000106 [-0.00112, 0.000912] | 1.00015 [0.99904, 1.00126]  |
| 150 MeV | equivalent | 1.00000 [0.99890, 1.00110] | 0.000425 [-0.00264, 0.00349] | 1.00028 [0.99899, 1.00158]  |

## B-picc against B-up (confirmatory)

- Joint claim (all 34 endpoints equivalent): **ESTABLISHED**
- Secondary joint claim (every endpoint other than rows where every run of both arms is zero is equivalent after Holm over all 34): **ESTABLISHED**; rows excluded: none
- Equivalent: 34 of 34 unadjusted, 34 after Holm; not established: 0

| endpoint | arm mean | reference mean | estimate | 90% interval | margin | outcome | Holm p |
|---|---|---|---|---|---|---|---|
| E100/lateral_39_5 | 4.2540 | 4.2612 | -0.00723 | [-0.0252, 0.0107] | [-0.5000, 0.5000] (difference) | equivalent | 8.61e-16 |
| E100/lateral_39_10 | 2.1171 | 2.1133 | 0.00374 | [-0.00972, 0.0172] | [-0.2000, 0.2000] (difference) | equivalent | 3.58e-12 |
| E100/lateral_39_20 | 1.1084 | 1.0994 | 0.00901 | [-0.00217, 0.0202] | [-0.2000, 0.2000] (difference) | equivalent | 5.18e-13 |
| E100/lateral_39_30 | 0.6840 | 0.6847 | -0.000765 | [-0.00611, 0.00458] | [-0.2000, 0.2000] (difference) | equivalent | 3.02e-17 |
| E100/lateral_39_50 | 0.2487 | 0.2479 | 0.000816 | [-0.00265, 0.00429] | [-0.2000, 0.2000] (difference) | equivalent | 5.55e-15 |
| E100/lateral_39_70 | 0.0834 | 0.0827 | 0.000773 | [-0.00184, 0.00338] | [-0.2000, 0.2000] (difference) | equivalent | 4.2e-20 |
| E100/lateral_61_5 | 5.9381 | 5.9403 | -0.00217 | [-0.0237, 0.0193] | [-0.5000, 0.5000] (difference) | equivalent | 1.3e-14 |
| E100/lateral_61_10 | 2.6570 | 2.6531 | 0.00393 | [-0.0103, 0.0181] | [-0.2000, 0.2000] (difference) | equivalent | 4.73e-12 |
| E100/lateral_61_20 | 1.1986 | 1.1944 | 0.00428 | [-0.005, 0.0136] | [-0.2000, 0.2000] (difference) | equivalent | 5.6e-14 |
| E100/lateral_61_30 | 0.7232 | 0.7104 | 0.0128 | [0.00526, 0.0203] | [-0.2000, 0.2000] (difference) | equivalent | 2.56e-14 |
| E100/lateral_61_50 | 0.2656 | 0.2618 | 0.00372 | [0.000232, 0.00721] | [-0.2000, 0.2000] (difference) | equivalent | 1.91e-18 |
| E100/lateral_61_70 | 0.0896 | 0.0879 | 0.00171 | [-0.00145, 0.00487] | [-0.2000, 0.2000] (difference) | equivalent | 1.32e-20 |
| E100/cax_39 | 50971.1470 | 50958.6147 | 1.0002 | [0.9993, 1.0012] | [0.9950, 1.0050] (ratio) | equivalent | 1.77e-06 |
| E100/cax_61 | 70270.4736 | 70252.2012 | 1.0003 | [0.9991, 1.0014] | [0.9950, 1.0050] (ratio) | equivalent | 9.34e-06 |
| E100/cax_max | 173690.9355 | 173549.5547 | 1.0008 | [0.9997, 1.0020] | [0.9950, 1.0050] (ratio) | equivalent | 2.9e-05 |
| E100/r80_mm | 75.7976 | 75.7976 | -4.89e-05 | [-0.000924, 0.000826] | [-0.5000, 0.5000] (difference) | equivalent | 4.97e-33 |
| E100/r20_mm | 78.3161 | 78.3159 | 0.000239 | [-0.000802, 0.00128] | [-0.5000, 0.5000] (difference) | equivalent | 1.52e-31 |
| E150/lateral_79_5 | 4.3589 | 4.3702 | -0.0113 | [-0.0342, 0.0116] | [-0.5000, 0.5000] (difference) | equivalent | 5.6e-14 |
| E150/lateral_79_10 | 2.6004 | 2.6228 | -0.0224 | [-0.0406, -0.00415] | [-0.2000, 0.2000] (difference) | equivalent | 2.92e-10 |
| E150/lateral_79_20 | 1.0434 | 1.0471 | -0.00377 | [-0.0164, 0.0089] | [-0.2000, 0.2000] (difference) | equivalent | 8.26e-12 |
| E150/lateral_79_30 | 0.4105 | 0.4059 | 0.00456 | [-0.00217, 0.0113] | [-0.2000, 0.2000] (difference) | equivalent | 5.18e-13 |
| E150/lateral_79_50 | 0.0559 | 0.0550 | 0.000869 | [-0.00165, 0.00339] | [-0.2000, 0.2000] (difference) | equivalent | 2.25e-21 |
| E150/lateral_79_70 | 0.00635 | 0.00624 | 0.000108 | [-0.00093, 0.00115] | [-0.2000, 0.2000] (difference) | equivalent | 5.83e-27 |
| E150/lateral_125_5 | 8.4214 | 8.4325 | -0.0111 | [-0.0356, 0.0134] | [-0.5000, 0.5000] (difference) | equivalent | 3.79e-12 |
| E150/lateral_125_10 | 3.6152 | 3.6261 | -0.0109 | [-0.0273, 0.00558] | [-0.2000, 0.2000] (difference) | equivalent | 1.59e-10 |
| E150/lateral_125_20 | 1.2813 | 1.2820 | -0.000693 | [-0.0101, 0.00876] | [-0.2000, 0.2000] (difference) | equivalent | 3.88e-12 |
| E150/lateral_125_30 | 0.4758 | 0.4791 | -0.00337 | [-0.0113, 0.00453] | [-0.2000, 0.2000] (difference) | equivalent | 1.27e-12 |
| E150/lateral_125_50 | 0.0581 | 0.0574 | 0.000695 | [-0.0016, 0.00299] | [-0.2000, 0.2000] (difference) | equivalent | 3.21e-21 |
| E150/lateral_125_70 | 0.00458 | 0.00452 | 5.65e-05 | [-0.000398, 0.000511] | [-0.2000, 0.2000] (difference) | equivalent | 2.04e-32 |
| E150/cax_79 | 36448.7949 | 36458.0591 | 0.9997 | [0.9984, 1.0011] | [0.9950, 1.0050] (ratio) | equivalent | 2.99e-05 |
| E150/cax_125 | 46420.0469 | 46437.5942 | 0.9996 | [0.9983, 1.0010] | [0.9950, 1.0050] (ratio) | equivalent | 2.99e-05 |
| E150/cax_max | 111974.7217 | 111967.8193 | 1.0001 | [0.9992, 1.0009] | [0.9950, 1.0050] (ratio) | equivalent | 3.35e-07 |
| E150/r80_mm | 156.1105 | 156.1132 | -0.00267 | [-0.00687, 0.00154] | [-0.5000, 0.5000] (difference) | equivalent | 1.22e-22 |
| E150/r20_mm | 159.6990 | 159.6978 | 0.00117 | [-0.00322, 0.00556] | [-0.5000, 0.5000] (difference) | equivalent | 3.96e-24 |

Descriptive only. Ratio of arm means (arm / reference) at 50 and 70 mm outside the edge, 95% interval:

| endpoint | arm mean | reference mean | ratio | 95% interval | note |
|---|---|---|---|---|---|
| E100/lateral_39_50 | 0.2487 | 0.2479 | 1.0033 | [0.9862, 1.0207] |  |
| E100/lateral_39_70 | 0.0834 | 0.0827 | 1.0094 | [0.9713, 1.0489] |  |
| E100/lateral_61_50 | 0.2656 | 0.2618 | 1.0142 | [0.9979, 1.0308] |  |
| E100/lateral_61_70 | 0.0896 | 0.0879 | 1.0195 | [0.9762, 1.0646] |  |
| E150/lateral_79_50 | 0.0559 | 0.0550 | 1.0158 | [0.9612, 1.0735] |  |
| E150/lateral_79_70 | 0.00635 | 0.00624 | 1.0173 | [0.8324, 1.2433] |  |
| E150/lateral_125_50 | 0.0581 | 0.0574 | 1.0121 | [0.9644, 1.0621] |  |
| E150/lateral_125_70 | 0.00458 | 0.00452 | 1.0125 | [0.8966, 1.1433] |  |

Descriptive only. Beside each maximum: the R80 difference and the ratio of the three-row mean around each run's maximum row. No cause is attributed.

| energy | maximum: outcome | maximum: ratio, 90% | R80 difference (mm), 90% | three-row mean: ratio, 95% |
|---|---|---|---|---|
| 100 MeV | equivalent | 1.00081 [0.99965, 1.00198] | -4.89e-05 [-0.000924, 0.000826] | 1.00078 [0.99950, 1.00205]  |
| 150 MeV | equivalent | 1.00006 [0.99921, 1.00091] | -0.00267 [-0.00687, 0.00154] | 0.99996 [0.99882, 1.00111]  |

## B-pgcc against B-picc (descriptive)

| endpoint | arm mean | reference mean | estimate | 90% interval | margin | outcome | Holm p |
|---|---|---|---|---|---|---|---|
| E100/lateral_39_5 | 4.2623 | 4.2540 | 0.00837 | [-0.0113, 0.0281] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E100/lateral_39_10 | 2.1178 | 2.1171 | 0.000761 | [-0.0122, 0.0137] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_39_20 | 1.1000 | 1.1084 | -0.00845 | [-0.0178, 0.000899] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_39_30 | 0.6814 | 0.6840 | -0.00257 | [-0.00831, 0.00316] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_39_50 | 0.2471 | 0.2487 | -0.0016 | [-0.00624, 0.00303] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_39_70 | 0.0852 | 0.0834 | 0.00178 | [-0.00118, 0.00473] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_61_5 | 5.9437 | 5.9381 | 0.00557 | [-0.0190, 0.0301] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E100/lateral_61_10 | 2.6589 | 2.6570 | 0.00195 | [-0.0125, 0.0164] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_61_20 | 1.1869 | 1.1986 | -0.0117 | [-0.0216, -0.00189] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_61_30 | 0.7126 | 0.7232 | -0.0105 | [-0.0166, -0.00451] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_61_50 | 0.2663 | 0.2656 | 0.000778 | [-0.00205, 0.00361] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/lateral_61_70 | 0.0893 | 0.0896 | -0.000314 | [-0.00338, 0.00275] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E100/cax_39 | 50966.1626 | 50971.1470 | 0.9999 | [0.9988, 1.0010] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E100/cax_61 | 70255.7920 | 70270.4736 | 0.9998 | [0.9989, 1.0007] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E100/cax_max | 173556.9570 | 173690.9355 | 0.9992 | [0.9982, 1.0003] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E100/r80_mm | 75.7975 | 75.7976 | -5.67e-05 | [-0.001, 0.000889] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E100/r20_mm | 78.3160 | 78.3161 | -0.000158 | [-0.00119, 0.000875] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E150/lateral_79_5 | 4.3548 | 4.3589 | -0.00412 | [-0.0223, 0.0140] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E150/lateral_79_10 | 2.6187 | 2.6004 | 0.0183 | [-0.000165, 0.0368] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_79_20 | 1.0475 | 1.0434 | 0.00416 | [-0.00513, 0.0135] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_79_30 | 0.4046 | 0.4105 | -0.00586 | [-0.0109, -0.00086] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_79_50 | 0.0538 | 0.0559 | -0.00205 | [-0.00517, 0.00108] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_79_70 | 0.00636 | 0.00635 | 1.33e-05 | [-0.000898, 0.000925] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_125_5 | 8.4438 | 8.4214 | 0.0224 | [-0.00771, 0.0525] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E150/lateral_125_10 | 3.6145 | 3.6152 | -0.000686 | [-0.0191, 0.0178] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_125_20 | 1.2870 | 1.2813 | 0.00574 | [-0.00171, 0.0132] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_125_30 | 0.4767 | 0.4758 | 0.00092 | [-0.00887, 0.0107] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_125_50 | 0.0572 | 0.0581 | -0.000918 | [-0.00308, 0.00125] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/lateral_125_70 | 0.00457 | 0.00458 | -6.28e-06 | [-0.000747, 0.000734] | [-0.2000, 0.2000] (difference) | descriptive |  |
| E150/cax_79 | 36444.1035 | 36448.7949 | 0.9999 | [0.9987, 1.0010] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E150/cax_125 | 46440.3403 | 46420.0469 | 1.0004 | [0.9989, 1.0020] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E150/cax_max | 111967.8330 | 111974.7217 | 0.9999 | [0.9989, 1.0010] | [0.9950, 1.0050] (ratio) | descriptive |  |
| E150/r80_mm | 156.1136 | 156.1105 | 0.00309 | [-0.000859, 0.00704] | [-0.5000, 0.5000] (difference) | descriptive |  |
| E150/r20_mm | 159.6991 | 159.6990 | 0.00012 | [-0.00327, 0.00351] | [-0.5000, 0.5000] (difference) | descriptive |  |

Descriptive only. Ratio of arm means (arm / reference) at 50 and 70 mm outside the edge, 95% interval:

| endpoint | arm mean | reference mean | ratio | 95% interval | note |
|---|---|---|---|---|---|
| E100/lateral_39_50 | 0.2471 | 0.2487 | 0.9936 | [0.9711, 1.0165] |  |
| E100/lateral_39_70 | 0.0852 | 0.0834 | 1.0213 | [0.9787, 1.0657] |  |
| E100/lateral_61_50 | 0.2663 | 0.2656 | 1.0029 | [0.9900, 1.0160] |  |
| E100/lateral_61_70 | 0.0893 | 0.0896 | 0.9965 | [0.9558, 1.0389] |  |
| E150/lateral_79_50 | 0.0538 | 0.0559 | 0.9634 | [0.8984, 1.0330] |  |
| E150/lateral_79_70 | 0.00636 | 0.00635 | 1.0021 | [0.8402, 1.1952] |  |
| E150/lateral_125_50 | 0.0572 | 0.0581 | 0.9842 | [0.9402, 1.0303] |  |
| E150/lateral_125_70 | 0.00457 | 0.00458 | 0.9986 | [0.8186, 1.2182] |  |

Descriptive only. Beside each maximum: the R80 difference and the ratio of the three-row mean around each run's maximum row. No cause is attributed.

| energy | maximum: outcome | maximum: ratio, 90% | R80 difference (mm), 90% | three-row mean: ratio, 95% |
|---|---|---|---|---|
| 100 MeV | descriptive | 0.99923 [0.99819, 1.00027] | -5.67e-05 [-0.001, 0.000889] | 0.99937 [0.99812, 1.00062]  |
| 150 MeV | descriptive | 0.99994 [0.99891, 1.00097] | 0.00309 [-0.000859, 0.00704] | 1.00032 [0.99907, 1.00157]  |

