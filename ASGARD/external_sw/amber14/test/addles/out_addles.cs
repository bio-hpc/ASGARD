add_les> ~ demo input for addles
add_les> ~
add_les> ~ addles input file for leucine dipeptide
add_les> ~ this makes a topology and coordinate file usable by SANDER
add_les> ~
add_les> ~ open input topology file
add_les> ~
add_les> file rprm name=(leu.dipep.prm) read
  The following unit number was assigned    26
|  INFO: Old style PARM file read

 Checking topology sizes against compiled limits
 Checking topology sizes against compiled limits
add_les> ~
add_les> ~ open input coordinates and velocities
add_les> ~
add_les> file rcvd name=(equ.crd) read
  The following unit number was assigned    27
 Coordinates and velocities from unit           27
for leu -> glu                                                                  
 Reading coordinates from input file
 Reading velocities from input file
add_les> ~
add_les> ~ open output topology
add_les> ~
add_les> file wprm name=(les.prm) wovr
  The following unit number was assigned    28
add_les> ~
add_les> ~ open output coordinates
add_les> ~
add_les> file wcrd name=(les.equ.crd) wovr
  The following unit number was assigned    29
add_les> ~
add_les> ~ all done with files, start processing commands
add_les> ~
add_les> action
add_les> ~
add_les> ~ first copy the leucine residue, 5 times
add_les> ~
add_les> spac numc=5 pick chem mono LEU done
 picking from           31  particles
atom, #, orig atom #, copies to make: N        7     7     5
atom, #, orig atom #, copies to make: H        8     8     5
atom, #, orig atom #, copies to make: CA       9     9     5
atom, #, orig atom #, copies to make: HA      10    10     5
atom, #, orig atom #, copies to make: CB      11    11     5
atom, #, orig atom #, copies to make: HB2     12    12     5
atom, #, orig atom #, copies to make: HB3     13    13     5
atom, #, orig atom #, copies to make: CG      14    14     5
atom, #, orig atom #, copies to make: HG      15    15     5
atom, #, orig atom #, copies to make: CD1     16    16     5
atom, #, orig atom #, copies to make: HD11    17    17     5
atom, #, orig atom #, copies to make: HD12    18    18     5
atom, #, orig atom #, copies to make: HD13    19    19     5
atom, #, orig atom #, copies to make: CD2     20    20     5
atom, #, orig atom #, copies to make: HD21    21    21     5
atom, #, orig atom #, copies to make: HD22    22    22     5
atom, #, orig atom #, copies to make: HD23    23    23     5
atom, #, orig atom #, copies to make: C       24    24     5
atom, #, orig atom #, copies to make: O       25    25     5
 Picked           19  particles
 Making            5  copies
Modifying velocities, factor is  2.37713
Modifying velocities, factor is  2.15761
Modifying velocities, factor is  2.10761
Modifying velocities, factor is  2.32795
Modifying velocities, factor is  2.33143
Modifying velocities, factor is  2.12593
Modifying velocities, factor is  2.15251
Modifying velocities, factor is  2.08598
Modifying velocities, factor is  2.28476
Modifying velocities, factor is  2.30802
Modifying velocities, factor is  2.38169
Modifying velocities, factor is  2.01381
Modifying velocities, factor is  2.14815
Modifying velocities, factor is  2.28375
Modifying velocities, factor is  2.11744
Modifying velocities, factor is  2.35496
Modifying velocities, factor is  2.35380
Modifying velocities, factor is  2.34510
Modifying velocities, factor is  2.08547
Modifying velocities, factor is  2.08056
Modifying velocities, factor is  2.33660
Modifying velocities, factor is  2.34147
Modifying velocities, factor is  2.17190
Modifying velocities, factor is  2.24333
Modifying velocities, factor is  2.41085
Modifying velocities, factor is  2.32327
Modifying velocities, factor is  2.23116
Modifying velocities, factor is  2.02266
Modifying velocities, factor is  2.22357
Modifying velocities, factor is  2.21988
Modifying velocities, factor is  2.02459
Modifying velocities, factor is  2.07497
Modifying velocities, factor is  2.05281
Modifying velocities, factor is  2.18657
Modifying velocities, factor is  2.18223
Modifying velocities, factor is  2.30448
Modifying velocities, factor is  2.17185
Modifying velocities, factor is  2.34720
Modifying velocities, factor is  2.08100
Modifying velocities, factor is  2.16038
Modifying velocities, factor is  2.30874
Modifying velocities, factor is  2.30560
Modifying velocities, factor is  2.21813
Modifying velocities, factor is  2.15837
Modifying velocities, factor is  2.15137
Modifying velocities, factor is  2.28781
Modifying velocities, factor is  2.13297
Modifying velocities, factor is  2.14312
Modifying velocities, factor is  2.06339
Modifying velocities, factor is  2.10416
Modifying velocities, factor is  2.26717
Modifying velocities, factor is  2.17820
Modifying velocities, factor is  2.36754
Modifying velocities, factor is  2.01279
Modifying velocities, factor is  2.34526
Modifying velocities, factor is  2.38478
Modifying velocities, factor is  2.17955
Modifying velocities, factor is  2.35033
Modifying velocities, factor is  2.14489
Modifying velocities, factor is  2.38189
Modifying velocities, factor is  2.10897
Modifying velocities, factor is  2.08195
Modifying velocities, factor is  2.07073
Modifying velocities, factor is  2.04585
Modifying velocities, factor is  2.36617
Modifying velocities, factor is  2.20236
Modifying velocities, factor is  2.43355
Modifying velocities, factor is  2.15602
Modifying velocities, factor is  2.41184
Modifying velocities, factor is  2.09164
Modifying velocities, factor is  2.28798
Modifying velocities, factor is  2.28406
Modifying velocities, factor is  2.20144
Modifying velocities, factor is  2.09794
Modifying velocities, factor is  2.33216
Modifying velocities, factor is  2.04318
Modifying velocities, factor is  2.09362
Modifying velocities, factor is  2.31949
Modifying velocities, factor is  2.21901
Modifying velocities, factor is  2.42687
Modifying velocities, factor is  2.35985
Modifying velocities, factor is  2.24695
Modifying velocities, factor is  2.25563
Modifying velocities, factor is  2.24649
Modifying velocities, factor is  2.22331
Modifying velocities, factor is  2.28254
Modifying velocities, factor is  2.24216
Modifying velocities, factor is  2.06708
Modifying velocities, factor is  2.30824
Modifying velocities, factor is  2.21472
Modifying velocities, factor is  2.09754
Modifying velocities, factor is  2.06947
Modifying velocities, factor is  2.15389
Modifying velocities, factor is  2.45212
Modifying velocities, factor is  2.30771
there were     31 particles; currently    107 particles
there were     18 nbonh bonds, now there are     62
there were     12 nbona bonds, now there are     48
there were      0 nbper bonds, now there are      0
there were     39 ntheth angles, now there are    143
there were     15 ntheta angles, now there are     71
there were      0 ngper angles, now there are      0
there were     56 nphih torsions, now there are    256
there were     16 nphia torsions, now there are     80
there were      0 ndper torsions, now there are      0
processing exclusion list 
finished creating LES subspace 
 Checking topology sizes against compiled limits
add_les> ~
add_les> ~ now  make an extra 2 copies of the atoms
add_les> ~ 14-23 (original numbering), just to check multi-level subspacing
add_les> ~
add_les> ~ this means each of the copies from above will have
add_les> ~ 2 copies made, resulting in 10 total copies
add_les> ~ of this section of the side chain
add_les> ~
add_les> spac numc=2 pick #prt 14 23 done
 picking from          107  particles
atom, #, orig atom #, copies to make: CG      42    14     2
atom, #, orig atom #, copies to make: CG      43    14     2
atom, #, orig atom #, copies to make: CG      44    14     2
atom, #, orig atom #, copies to make: CG      45    14     2
atom, #, orig atom #, copies to make: CG      46    14     2
atom, #, orig atom #, copies to make: HG      47    15     2
atom, #, orig atom #, copies to make: HG      48    15     2
atom, #, orig atom #, copies to make: HG      49    15     2
atom, #, orig atom #, copies to make: HG      50    15     2
atom, #, orig atom #, copies to make: HG      51    15     2
atom, #, orig atom #, copies to make: CD1     52    16     2
atom, #, orig atom #, copies to make: CD1     53    16     2
atom, #, orig atom #, copies to make: CD1     54    16     2
atom, #, orig atom #, copies to make: CD1     55    16     2
atom, #, orig atom #, copies to make: CD1     56    16     2
atom, #, orig atom #, copies to make: HD11    57    17     2
atom, #, orig atom #, copies to make: HD11    58    17     2
atom, #, orig atom #, copies to make: HD11    59    17     2
atom, #, orig atom #, copies to make: HD11    60    17     2
atom, #, orig atom #, copies to make: HD11    61    17     2
atom, #, orig atom #, copies to make: HD12    62    18     2
atom, #, orig atom #, copies to make: HD12    63    18     2
atom, #, orig atom #, copies to make: HD12    64    18     2
atom, #, orig atom #, copies to make: HD12    65    18     2
atom, #, orig atom #, copies to make: HD12    66    18     2
atom, #, orig atom #, copies to make: HD13    67    19     2
atom, #, orig atom #, copies to make: HD13    68    19     2
atom, #, orig atom #, copies to make: HD13    69    19     2
atom, #, orig atom #, copies to make: HD13    70    19     2
atom, #, orig atom #, copies to make: HD13    71    19     2
atom, #, orig atom #, copies to make: CD2     72    20     2
atom, #, orig atom #, copies to make: CD2     73    20     2
atom, #, orig atom #, copies to make: CD2     74    20     2
atom, #, orig atom #, copies to make: CD2     75    20     2
atom, #, orig atom #, copies to make: CD2     76    20     2
atom, #, orig atom #, copies to make: HD21    77    21     2
atom, #, orig atom #, copies to make: HD21    78    21     2
atom, #, orig atom #, copies to make: HD21    79    21     2
atom, #, orig atom #, copies to make: HD21    80    21     2
atom, #, orig atom #, copies to make: HD21    81    21     2
atom, #, orig atom #, copies to make: HD22    82    22     2
atom, #, orig atom #, copies to make: HD22    83    22     2
atom, #, orig atom #, copies to make: HD22    84    22     2
atom, #, orig atom #, copies to make: HD22    85    22     2
atom, #, orig atom #, copies to make: HD22    86    22     2
atom, #, orig atom #, copies to make: HD23    87    23     2
atom, #, orig atom #, copies to make: HD23    88    23     2
atom, #, orig atom #, copies to make: HD23    89    23     2
atom, #, orig atom #, copies to make: HD23    90    23     2
atom, #, orig atom #, copies to make: HD23    91    23     2
 Picked           50  particles
 Making            2  copies
Modifying velocities, factor is  1.52331
Modifying velocities, factor is  1.45082
Modifying velocities, factor is  1.53093
Modifying velocities, factor is  1.34190
Modifying velocities, factor is  1.29311
Modifying velocities, factor is  1.45418
Modifying velocities, factor is  1.42363
Modifying velocities, factor is  1.34237
Modifying velocities, factor is  1.36405
Modifying velocities, factor is  1.45349
Modifying velocities, factor is  1.47767
Modifying velocities, factor is  1.40403
Modifying velocities, factor is  1.43154
Modifying velocities, factor is  1.51068
Modifying velocities, factor is  1.40995
Modifying velocities, factor is  1.48189
Modifying velocities, factor is  1.28646
Modifying velocities, factor is  1.43426
Modifying velocities, factor is  1.32211
Modifying velocities, factor is  1.51705
Modifying velocities, factor is  1.47640
Modifying velocities, factor is  1.32528
Modifying velocities, factor is  1.34807
Modifying velocities, factor is  1.44273
Modifying velocities, factor is  1.49007
Modifying velocities, factor is  1.49218
Modifying velocities, factor is  1.29180
Modifying velocities, factor is  1.42460
Modifying velocities, factor is  1.28666
Modifying velocities, factor is  1.30158
Modifying velocities, factor is  1.50430
Modifying velocities, factor is  1.44239
Modifying velocities, factor is  1.53926
Modifying velocities, factor is  1.35353
Modifying velocities, factor is  1.54106
Modifying velocities, factor is  1.39173
Modifying velocities, factor is  1.42422
Modifying velocities, factor is  1.39653
Modifying velocities, factor is  1.28061
Modifying velocities, factor is  1.55111
Modifying velocities, factor is  1.31021
Modifying velocities, factor is  1.46778
Modifying velocities, factor is  1.31836
Modifying velocities, factor is  1.42125
Modifying velocities, factor is  1.28576
Modifying velocities, factor is  1.34950
Modifying velocities, factor is  1.39497
Modifying velocities, factor is  1.42820
Modifying velocities, factor is  1.38736
Modifying velocities, factor is  1.43542
Modifying velocities, factor is  1.36623
Modifying velocities, factor is  1.35037
Modifying velocities, factor is  1.38757
Modifying velocities, factor is  1.50099
Modifying velocities, factor is  1.33522
Modifying velocities, factor is  1.49995
Modifying velocities, factor is  1.46176
Modifying velocities, factor is  1.28561
Modifying velocities, factor is  1.28275
Modifying velocities, factor is  1.43733
Modifying velocities, factor is  1.30365
Modifying velocities, factor is  1.43768
Modifying velocities, factor is  1.37919
Modifying velocities, factor is  1.28842
Modifying velocities, factor is  1.47227
Modifying velocities, factor is  1.48872
Modifying velocities, factor is  1.44040
Modifying velocities, factor is  1.27831
Modifying velocities, factor is  1.42004
Modifying velocities, factor is  1.34108
Modifying velocities, factor is  1.37676
Modifying velocities, factor is  1.41921
Modifying velocities, factor is  1.40193
Modifying velocities, factor is  1.48655
Modifying velocities, factor is  1.40577
Modifying velocities, factor is  1.33593
Modifying velocities, factor is  1.51022
Modifying velocities, factor is  1.30913
Modifying velocities, factor is  1.40634
Modifying velocities, factor is  1.34992
Modifying velocities, factor is  1.38219
Modifying velocities, factor is  1.40068
Modifying velocities, factor is  1.43928
Modifying velocities, factor is  1.28327
Modifying velocities, factor is  1.43366
Modifying velocities, factor is  1.51972
Modifying velocities, factor is  1.52069
Modifying velocities, factor is  1.28728
Modifying velocities, factor is  1.50888
Modifying velocities, factor is  1.28255
Modifying velocities, factor is  1.45742
Modifying velocities, factor is  1.53017
Modifying velocities, factor is  1.43052
Modifying velocities, factor is  1.39266
Modifying velocities, factor is  1.44105
Modifying velocities, factor is  1.53436
Modifying velocities, factor is  1.38000
Modifying velocities, factor is  1.38897
Modifying velocities, factor is  1.42907
Modifying velocities, factor is  1.40372
there were    107 particles; currently    157 particles
there were     62 nbonh bonds, now there are     97
there were     48 nbona bonds, now there are     63
there were      0 nbper bonds, now there are      0
there were    143 ntheth angles, now there are    228
there were     71 ntheta angles, now there are     91
there were      0 ngper angles, now there are      0
there were    256 nphih torsions, now there are    386
there were     80 nphia torsions, now there are    100
there were      0 ndper torsions, now there are      0
processing exclusion list 
finished creating LES subspace 
 Checking topology sizes against compiled limits
add_les> ~
add_les> ~ this end line needs to be present!
add_les> ~
add_les> *EOD
 Finished reading subspace definitions. 
 Looking for unique atom and covalent types
 There were            8  bond types, now there are           18
 There were           20  angle types, now there are           34
 multi-term torsion          24           6
 multi-term torsion          16          19
 multi-term torsion          16          20
 multi-term torsion          81          27
 multi-term torsion          81          28
 multi-term torsion          86          33
 There were           19  dihedral types, now there are           35
 There were            7  atom typesNow there are           16
 MAX # ATOMS IN A SINGLE RESIDUE =          145
  ATOM  origpt
ATOM      1 original atom #      1
ATOM      2 original atom #      2
ATOM      3 original atom #      3
ATOM      4 original atom #      4
ATOM      5 original atom #      5
ATOM      6 original atom #      6
ATOM      7 original atom #      7
ATOM      8 original atom #      7
ATOM      9 original atom #      7
ATOM     10 original atom #      7
ATOM     11 original atom #      7
ATOM     12 original atom #      8
ATOM     13 original atom #      8
ATOM     14 original atom #      8
ATOM     15 original atom #      8
ATOM     16 original atom #      8
ATOM     17 original atom #      9
ATOM     18 original atom #      9
ATOM     19 original atom #      9
ATOM     20 original atom #      9
ATOM     21 original atom #      9
ATOM     22 original atom #     10
ATOM     23 original atom #     10
ATOM     24 original atom #     10
ATOM     25 original atom #     10
ATOM     26 original atom #     10
ATOM     27 original atom #     11
ATOM     28 original atom #     11
ATOM     29 original atom #     11
ATOM     30 original atom #     11
ATOM     31 original atom #     11
ATOM     32 original atom #     12
ATOM     33 original atom #     12
ATOM     34 original atom #     12
ATOM     35 original atom #     12
ATOM     36 original atom #     12
ATOM     37 original atom #     13
ATOM     38 original atom #     13
ATOM     39 original atom #     13
ATOM     40 original atom #     13
ATOM     41 original atom #     13
ATOM     42 original atom #     14
ATOM     43 original atom #     14
ATOM     44 original atom #     14
ATOM     45 original atom #     14
ATOM     46 original atom #     14
ATOM     47 original atom #     14
ATOM     48 original atom #     14
ATOM     49 original atom #     14
ATOM     50 original atom #     14
ATOM     51 original atom #     14
ATOM     52 original atom #     15
ATOM     53 original atom #     15
ATOM     54 original atom #     15
ATOM     55 original atom #     15
ATOM     56 original atom #     15
ATOM     57 original atom #     15
ATOM     58 original atom #     15
ATOM     59 original atom #     15
ATOM     60 original atom #     15
ATOM     61 original atom #     15
ATOM     62 original atom #     16
ATOM     63 original atom #     16
ATOM     64 original atom #     16
ATOM     65 original atom #     16
ATOM     66 original atom #     16
ATOM     67 original atom #     16
ATOM     68 original atom #     16
ATOM     69 original atom #     16
ATOM     70 original atom #     16
ATOM     71 original atom #     16
ATOM     72 original atom #     17
ATOM     73 original atom #     17
ATOM     74 original atom #     17
ATOM     75 original atom #     17
ATOM     76 original atom #     17
ATOM     77 original atom #     17
ATOM     78 original atom #     17
ATOM     79 original atom #     17
ATOM     80 original atom #     17
ATOM     81 original atom #     17
ATOM     82 original atom #     18
ATOM     83 original atom #     18
ATOM     84 original atom #     18
ATOM     85 original atom #     18
ATOM     86 original atom #     18
ATOM     87 original atom #     18
ATOM     88 original atom #     18
ATOM     89 original atom #     18
ATOM     90 original atom #     18
ATOM     91 original atom #     18
ATOM     92 original atom #     19
ATOM     93 original atom #     19
ATOM     94 original atom #     19
ATOM     95 original atom #     19
ATOM     96 original atom #     19
ATOM     97 original atom #     19
ATOM     98 original atom #     19
ATOM     99 original atom #     19
ATOM    100 original atom #     19
ATOM    101 original atom #     19
ATOM    102 original atom #     20
ATOM    103 original atom #     20
ATOM    104 original atom #     20
ATOM    105 original atom #     20
ATOM    106 original atom #     20
ATOM    107 original atom #     20
ATOM    108 original atom #     20
ATOM    109 original atom #     20
ATOM    110 original atom #     20
ATOM    111 original atom #     20
ATOM    112 original atom #     21
ATOM    113 original atom #     21
ATOM    114 original atom #     21
ATOM    115 original atom #     21
ATOM    116 original atom #     21
ATOM    117 original atom #     21
ATOM    118 original atom #     21
ATOM    119 original atom #     21
ATOM    120 original atom #     21
ATOM    121 original atom #     21
ATOM    122 original atom #     22
ATOM    123 original atom #     22
ATOM    124 original atom #     22
ATOM    125 original atom #     22
ATOM    126 original atom #     22
ATOM    127 original atom #     22
ATOM    128 original atom #     22
ATOM    129 original atom #     22
ATOM    130 original atom #     22
ATOM    131 original atom #     22
ATOM    132 original atom #     23
ATOM    133 original atom #     23
ATOM    134 original atom #     23
ATOM    135 original atom #     23
ATOM    136 original atom #     23
ATOM    137 original atom #     23
ATOM    138 original atom #     23
ATOM    139 original atom #     23
ATOM    140 original atom #     23
ATOM    141 original atom #     23
ATOM    142 original atom #     24
ATOM    143 original atom #     24
ATOM    144 original atom #     24
ATOM    145 original atom #     24
ATOM    146 original atom #     24
ATOM    147 original atom #     25
ATOM    148 original atom #     25
ATOM    149 original atom #     25
ATOM    150 original atom #     25
ATOM    151 original atom #     25
ATOM    152 original atom #     26
ATOM    153 original atom #     27
ATOM    154 original atom #     28
ATOM    155 original atom #     29
ATOM    156 original atom #     30
ATOM    157 original atom #     31
 WARNING: atom           42
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           43
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           44
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           45
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           46
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           47
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           48
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           49
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           50
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           51
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           52
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           53
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           54
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           55
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           56
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           57
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           58
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           59
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           60
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           61
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           62
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           63
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           64
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           65
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           66
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           67
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           68
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           69
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           70
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           71
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           72
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           73
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           74
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           75
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           76
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           77
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           78
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           79
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           80
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           81
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           82
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           83
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           84
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           85
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           86
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           87
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           88
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           89
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           90
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           91
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           92
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           93
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           94
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           95
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           96
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           97
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           98
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom           99
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          100
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          101
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          102
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          103
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          104
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          105
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          106
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          107
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          108
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          109
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          110
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          111
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          112
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          113
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          114
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          115
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          116
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          117
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          118
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          119
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          120
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          121
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          122
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          123
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          124
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          125
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          126
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          127
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          128
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          129
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          130
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          131
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          132
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          133
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          134
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          135
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          136
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          137
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          138
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          139
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          140
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 WARNING: atom          141
 is in more than 1 LES subspace!
 NMR restraints may not be supported.
 Writing coordinates to output file
 Writing velocities to output file
 Successful completion of ADDLES
