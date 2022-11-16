! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  UFF VDW parameters (dependent only on atomic number not atom type)     |#
! #|                                                                         |#
! #|  [x] = Angstroms   [d] = kcal/mol                                       |#
! #|                                                                         |#
! #|  A. K. Rappe, C. J. Casewit, K. S. Colwell, W. A. Goddard III, and      |#
! #|  W. M. Skiff, "UFF, a Full Periodic Table Force Field for Molecular     |#
! #|  Mechanics and Molecular Dynamics Simulations", J. Am. Chem. Soc., 114  |#
! #|  (1992), 10024-10035.                                                   |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   integer, parameter :: UFFvdw_size = 103
   type vdw_UFF
      _REAL_ :: x
      _REAL_ :: d
      _REAL_ :: z
   end type vdw_UFF

   type( vdw_UFF ) UFFvdw(UFFvdw_size)


   UFFvdw( 1)%x = 2.886; UFFvdw( 1)%d = 0.044; UFFvdw( 1)%z = 12.000
   UFFvdw( 2)%x = 2.362; UFFvdw( 2)%d = 0.056; UFFvdw( 2)%z = 15.240
   UFFvdw( 3)%x = 2.451; UFFvdw( 3)%d = 0.025; UFFvdw( 3)%z = 12.000
   UFFvdw( 4)%x = 2.745; UFFvdw( 4)%d = 0.085; UFFvdw( 4)%z = 12.000
   UFFvdw( 5)%x = 4.083; UFFvdw( 5)%d = 0.180; UFFvdw( 5)%z = 12.052 
   UFFvdw( 6)%x = 3.851; UFFvdw( 6)%d = 0.105; UFFvdw( 6)%z = 12.730
   UFFvdw( 7)%x = 3.660; UFFvdw( 7)%d = 0.069; UFFvdw( 7)%z = 13.407
   UFFvdw( 8)%x = 3.500; UFFvdw( 8)%d = 0.060; UFFvdw( 8)%z = 14.085
   UFFvdw( 9)%x = 3.364; UFFvdw( 9)%d = 0.050; UFFvdw( 9)%z = 14.762
   UFFvdw(10)%x = 3.243; UFFvdw(10)%d = 0.042; UFFvdw(10)%z = 15.440
   UFFvdw(11)%x = 2.983; UFFvdw(11)%d = 0.030; UFFvdw(11)%z = 12.000
   UFFvdw(12)%x = 3.021; UFFvdw(12)%d = 0.111; UFFvdw(12)%z = 12.000
   UFFvdw(13)%x = 4.499; UFFvdw(13)%d = 0.505; UFFvdw(13)%z = 11.278
   UFFvdw(14)%x = 4.295; UFFvdw(14)%d = 0.402; UFFvdw(14)%z = 12.175
   UFFvdw(15)%x = 4.147; UFFvdw(15)%d = 0.305; UFFvdw(15)%z = 13.072
   UFFvdw(16)%x = 4.035; UFFvdw(16)%d = 0.274; UFFvdw(16)%z = 13.969
   UFFvdw(17)%x = 3.947; UFFvdw(17)%d = 0.227; UFFvdw(17)%z = 14.866
   UFFvdw(18)%x = 3.868; UFFvdw(18)%d = 0.185; UFFvdw(18)%z = 15.763
   UFFvdw(19)%x = 3.812; UFFvdw(19)%d = 0.035; UFFvdw(19)%z = 12.000
   UFFvdw(20)%x = 3.399; UFFvdw(20)%d = 0.238; UFFvdw(20)%z = 12.000
   UFFvdw(21)%x = 3.295; UFFvdw(21)%d = 0.019; UFFvdw(21)%z = 12.000
   UFFvdw(22)%x = 3.175; UFFvdw(22)%d = 0.017; UFFvdw(22)%z = 12.000
   UFFvdw(23)%x = 3.144; UFFvdw(23)%d = 0.016; UFFvdw(23)%z = 12.000
   UFFvdw(24)%x = 3.023; UFFvdw(24)%d = 0.015; UFFvdw(24)%z = 12.000
   UFFvdw(25)%x = 2.961; UFFvdw(25)%d = 0.013; UFFvdw(25)%z = 12.000

   UFFvdw(26)%x = 2.912; UFFvdw(26)%d = 0.013; UFFvdw(26)%z = 12.000
   UFFvdw(27)%x = 2.872; UFFvdw(27)%d = 0.014; UFFvdw(27)%z = 12.000
   UFFvdw(28)%x = 2.834; UFFvdw(28)%d = 0.015; UFFvdw(28)%z = 12.000
   UFFvdw(29)%x = 3.495; UFFvdw(29)%d = 0.005; UFFvdw(29)%z = 12.000
   UFFvdw(30)%x = 2.763; UFFvdw(30)%d = 0.124; UFFvdw(30)%z = 12.000
   UFFvdw(31)%x = 4.383; UFFvdw(31)%d = 0.415; UFFvdw(31)%z = 11.000
   UFFvdw(32)%x = 4.280; UFFvdw(32)%d = 0.379; UFFvdw(32)%z = 12.000
   UFFvdw(33)%x = 4.230; UFFvdw(33)%d = 0.309; UFFvdw(33)%z = 13.000
   UFFvdw(34)%x = 4.205; UFFvdw(34)%d = 0.291; UFFvdw(34)%z = 14.000
   UFFvdw(35)%x = 4.189; UFFvdw(35)%d = 0.251; UFFvdw(35)%z = 15.000
   UFFvdw(36)%x = 4.141; UFFvdw(36)%d = 0.220; UFFvdw(36)%z = 16.000
   UFFvdw(37)%x = 4.114; UFFvdw(37)%d = 0.040; UFFvdw(37)%z = 12.000
   UFFvdw(38)%x = 3.641; UFFvdw(38)%d = 0.235; UFFvdw(38)%z = 12.000
   UFFvdw(39)%x = 3.345; UFFvdw(39)%d = 0.072; UFFvdw(39)%z = 12.000
   UFFvdw(40)%x = 3.124; UFFvdw(40)%d = 0.069; UFFvdw(40)%z = 12.000
   UFFvdw(41)%x = 3.165; UFFvdw(41)%d = 0.059; UFFvdw(41)%z = 12.000
   UFFvdw(42)%x = 3.052; UFFvdw(42)%d = 0.056; UFFvdw(42)%z = 12.000
   UFFvdw(43)%x = 2.998; UFFvdw(43)%d = 0.048; UFFvdw(43)%z = 12.000
   UFFvdw(44)%x = 2.963; UFFvdw(44)%d = 0.056; UFFvdw(44)%z = 12.000
   UFFvdw(45)%x = 2.929; UFFvdw(45)%d = 0.053; UFFvdw(45)%z = 12.000
   UFFvdw(46)%x = 2.899; UFFvdw(46)%d = 0.048; UFFvdw(46)%z = 12.000
   UFFvdw(47)%x = 3.148; UFFvdw(47)%d = 0.036; UFFvdw(47)%z = 12.000
   UFFvdw(48)%x = 2.848; UFFvdw(48)%d = 0.228; UFFvdw(48)%z = 12.000
   UFFvdw(49)%x = 4.463; UFFvdw(49)%d = 0.599; UFFvdw(49)%z = 11.000
   UFFvdw(50)%x = 4.392; UFFvdw(50)%d = 0.567; UFFvdw(50)%z = 12.000

   UFFvdw(51)%x = 4.420; UFFvdw(51)%d = 0.449; UFFvdw(51)%z = 13.000
   UFFvdw(52)%x = 4.470; UFFvdw(52)%d = 0.398; UFFvdw(52)%z = 14.000
   UFFvdw(53)%x = 4.500; UFFvdw(53)%d = 0.339; UFFvdw(53)%z = 15.000
   UFFvdw(54)%x = 4.404; UFFvdw(54)%d = 0.332; UFFvdw(54)%z = 12.000
   UFFvdw(55)%x = 4.517; UFFvdw(55)%d = 0.045; UFFvdw(55)%z = 12.000
   UFFvdw(56)%x = 3.703; UFFvdw(56)%d = 0.364; UFFvdw(56)%z = 12.000
   UFFvdw(57)%x = 3.522; UFFvdw(57)%d = 0.017; UFFvdw(57)%z = 12.000
   UFFvdw(58)%x = 3.556; UFFvdw(58)%d = 0.130; UFFvdw(58)%z = 12.000
   UFFvdw(59)%x = 3.606; UFFvdw(59)%d = 0.010; UFFvdw(59)%z = 12.000
   UFFvdw(60)%x = 3.575; UFFvdw(60)%d = 0.010; UFFvdw(60)%z = 12.000
   UFFvdw(61)%x = 3.547; UFFvdw(61)%d = 0.009; UFFvdw(61)%z = 12.000
   UFFvdw(62)%x = 3.520; UFFvdw(62)%d = 0.008; UFFvdw(62)%z = 12.000
   UFFvdw(63)%x = 3.493; UFFvdw(63)%d = 0.008; UFFvdw(63)%z = 12.000
   UFFvdw(64)%x = 3.368; UFFvdw(64)%d = 0.009; UFFvdw(64)%z = 12.000
   UFFvdw(65)%x = 3.451; UFFvdw(65)%d = 0.007; UFFvdw(65)%z = 12.000
   UFFvdw(66)%x = 3.428; UFFvdw(66)%d = 0.007; UFFvdw(66)%z = 12.000
   UFFvdw(67)%x = 3.409; UFFvdw(67)%d = 0.007; UFFvdw(67)%z = 12.000
   UFFvdw(68)%x = 3.391; UFFvdw(68)%d = 0.007; UFFvdw(68)%z = 12.000
   UFFvdw(69)%x = 3.374; UFFvdw(69)%d = 0.006; UFFvdw(69)%z = 12.000
   UFFvdw(70)%x = 3.355; UFFvdw(70)%d = 0.228; UFFvdw(70)%z = 12.000
   UFFvdw(71)%x = 3.640; UFFvdw(71)%d = 0.041; UFFvdw(71)%z = 12.000
   UFFvdw(72)%x = 3.141; UFFvdw(72)%d = 0.072; UFFvdw(72)%z = 12.000
   UFFvdw(73)%x = 3.170; UFFvdw(73)%d = 0.081; UFFvdw(73)%z = 12.000
   UFFvdw(74)%x = 3.069; UFFvdw(74)%d = 0.067; UFFvdw(74)%z = 12.000
   UFFvdw(75)%x = 2.954; UFFvdw(75)%d = 0.066; UFFvdw(75)%z = 12.000

   UFFvdw(76)%x = 3.120; UFFvdw(76)%d = 0.037; UFFvdw(76)%z = 12.000
   UFFvdw(77)%x = 2.840; UFFvdw(77)%d = 0.073; UFFvdw(77)%z = 12.000
   UFFvdw(78)%x = 2.754; UFFvdw(78)%d = 0.080; UFFvdw(78)%z = 12.000
   UFFvdw(79)%x = 3.293; UFFvdw(79)%d = 0.039; UFFvdw(79)%z = 12.000
   UFFvdw(80)%x = 2.705; UFFvdw(80)%d = 0.385; UFFvdw(80)%z = 12.000
   UFFvdw(81)%x = 4.347; UFFvdw(81)%d = 0.680; UFFvdw(81)%z = 11.000
   UFFvdw(82)%x = 4.297; UFFvdw(82)%d = 0.663; UFFvdw(82)%z = 12.000
   UFFvdw(83)%x = 4.370; UFFvdw(83)%d = 0.518; UFFvdw(83)%z = 13.000
   UFFvdw(84)%x = 4.709; UFFvdw(84)%d = 0.325; UFFvdw(84)%z = 14.000
   UFFvdw(85)%x = 4.750; UFFvdw(85)%d = 0.284; UFFvdw(85)%z = 15.000
   UFFvdw(86)%x = 4.765; UFFvdw(86)%d = 0.248; UFFvdw(86)%z = 16.000
   UFFvdw(87)%x = 4.900; UFFvdw(87)%d = 0.050; UFFvdw(87)%z = 12.000
   UFFvdw(88)%x = 3.677; UFFvdw(88)%d = 0.404; UFFvdw(88)%z = 12.000
   UFFvdw(89)%x = 3.478; UFFvdw(89)%d = 0.033; UFFvdw(89)%z = 12.000
   UFFvdw(90)%x = 3.396; UFFvdw(90)%d = 0.026; UFFvdw(90)%z = 12.000
   UFFvdw(91)%x = 3.424; UFFvdw(91)%d = 0.022; UFFvdw(91)%z = 12.000
   UFFvdw(92)%x = 3.395; UFFvdw(92)%d = 0.022; UFFvdw(92)%z = 12.000
   UFFvdw(93)%x = 3.424; UFFvdw(93)%d = 0.019; UFFvdw(93)%z = 12.000
   UFFvdw(94)%x = 3.424; UFFvdw(94)%d = 0.016; UFFvdw(94)%z = 12.000
   UFFvdw(95)%x = 3.381; UFFvdw(95)%d = 0.014; UFFvdw(95)%z = 12.000
   UFFvdw(96)%x = 3.326; UFFvdw(96)%d = 0.013; UFFvdw(96)%z = 12.000
   UFFvdw(97)%x = 3.339; UFFvdw(97)%d = 0.013; UFFvdw(97)%z = 12.000
   UFFvdw(98)%x = 3.313; UFFvdw(98)%d = 0.013; UFFvdw(98)%z = 12.000
   UFFvdw(99)%x = 3.299; UFFvdw(99)%d = 0.012; UFFvdw(99)%z = 12.000
   UFFvdw(100)%x = 3.286; UFFvdw(100)%d = 0.012; UFFvdw(100)%z = 12.000

   UFFvdw(101)%x = 3.274; UFFvdw(101)%d = 0.011; UFFvdw(101)%z = 12.000
   UFFvdw(102)%x = 3.248; UFFvdw(102)%d = 0.011; UFFvdw(102)%z = 12.000
   UFFvdw(103)%x = 3.236; UFFvdw(103)%d = 0.011; UFFvdw(103)%z = 12.000

!  +---------------------------------------------------------------------------+
!  |  Convert [x] to Bohrs & [d] to Hartrees                                   |
!  +---------------------------------------------------------------------------+

!  UFFvdw(:)%x = UFFvdw(:)%x * A_TO_BOHRS
!  UFFvdw(:)%d = UFFvdw(:)%d / AU_TO_KCAL

   UFFvdw(:)%x = UFFvdw(:)%x / Bohr2Angstrom
   UFFvdw(:)%d = UFFvdw(:)%d / Hartree2Kcalmol


!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x


