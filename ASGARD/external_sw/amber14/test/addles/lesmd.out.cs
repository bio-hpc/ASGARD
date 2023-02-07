
          -------------------------------------------------------
          Amber 5  SANDER                   Scripps/UCSF 1997
          -------------------------------------------------------

|                                        Mon Aug 11 18:22:32 1997

  [-O]verwriting output

File Assignments:
|MDIN : lesmd.in                                                              
|MDOUT: lesmd.out                                                             
|INPCR: les.equ.crd                                                           
|PARM : les.prm                                                               
|RESTR: lesmd.restrt                                                          
|REFC : refc                                                                  
|MDVEL: lesmd.vel                                                             
|MDEN : lesmd.en                                                              
|MDCRD: lesmd.crd                                                             
|MDINF: lesmd.info                                                            


 Here is the input file:

leucine LES test case                                                          
 &cntrl                                                                        
                                                                               
  ntx    = 5, ntwx   = 10,   nsnb   = 10, idiel  = 1,    dielc = 1.0,          
  nrun   = 1, nstlim = 1000, ntpr   = 10, temp0  = 298., dt    = 0.001,        
  ntt    = 1, init   = 4,    ntc    = 1,  tol    = 0.0005,                     
  scee   = 1.2,                                                                
                                                                               
 /                                                                          
-------------------------------------------------------------------------------

leucine LES test case                                                           

| Reading &cntrl namelist w/ portable lib



   1.  RESOURCE   USE: 

 NATOM  =     157 NTYPES =      16 NBONH =      97 MBONA  =      63
 NTHETH =     228 MTHETA =      91 NPHIH =     386 MPHIA  =     100
 NHPARM =       0 NPARM  =       1 NNB   =    9856 NRES   =       3
 NBONA  =      63 NTHETA =      91 NPHIA =     100 NUMBND =      16
 NUMANG =      23 NPTRA  =      20 NATYP =       7 NPHB   =       0
 IFBOX  =       0 NMXRS  =     145 IFCAP =       0


|     Memory type      Allocated 
|     Real              108454
|     Hollerith           1261
|     Integer            96195 (static)

|     Max Nonbonded Pairs:   12325 packed  2 to a machine word
 LES parameters were found


   2.  CONTROL  DATA  FOR  THE  RUN

for leu -> glu                                                                  

     TIMLIM=  999999.   IREST =    0       IBELLY=    0
     KFORM =    1       ICHDNA=    0       IMIN  =    0
     IPOL  =    0       IEWALD=    0

     NTX   =    5       NTXO  =    1
     IG    =    71277   TEMPI =     0.00   HEAT  =    0.000

     NTB   =    0       IFTRES=    1       BOXX  =    0.000
     BOXY  =    0.000   BOXZ  =    0.000

     NRUN  =    1       NTT   =    1       TEMP0 =  298.000
     DTEMP =    0.000   TAUTP =    0.200   TAUTS =    0.200
     ISOLVP=    0       VLIMIT=    0.000

     NTP   =    0       PRES0 =    1.000   COMP  =   44.600
     TAUP  =    0.200   NPSCAL=    0

     NTCM  =    0       NSCM  = 9999999

     NSTLIM= 1000       INIT  =    4       NTU   =    1
     T     =    0.000   DT    =   0.00100

     NTC   =    1       TOL   =   0.00050  JFASTW =    0

     NTF   =    1       NTID  =    0       NTNB  =    1
     NSNB  =   10       IDIEL =    1       IMGSLT=    0
     IPRR  =    0       IPRW  =    0       ITRSLU=    1

     CUT   =    8.000   SCNB  =    2.000
     SCEE  =    1.200   DIELC =    1.000
     CUT2ND=   0.00000

     NTPR  =      10    NTWR  =      50    NTWX  =      10
     NTWV  =       0    NTWE  =       0    NTWXM =  999999
     NTWVM =  999999    NTWEM =  999999    IOUTFM=       0
     NTWPRT=       0    NTWPR0=       0

     NTR   =    0       NTRX  =    1
     TAUR  =   0.00000     NMROPT=    0     ISFTRP=    0
     RWELL =   1.00000     PENCUT=   0.10000

     IVCAP =    0       MATCAP=    0       FCAP  =    1.500

   OTHER DATA:

     IFCAP =    0       NATCAP=    0       CUTCAP=    0.000
     XCAP  =    0.000   YCAP  =    0.000   ZCAP  =    0.000

     NATOM =     157  NRES =      3

     Water definition for fast triangulated model:
     Resname = WAT ; Oxygen_name = O   ; Hyd1_name = H1  ; Hyd2_name = H2  

   3.  ATOMIC COORDINATES AND VELOCITIES

for leu -> glu                                                                  
 begin time read from input coords =     0.000 ps

 Number of triangulated 3-point waters found:        0

 Solute/solvent pointers:
     IPTSOL=    3       NATRCM=  157
     IPTRES=    0       IPTATM=    0
     NSPSOL=    0       NSPSTR=    0
     NSOLUT=  157       NATOM =  157

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =     1  TIME(PS) =    0.001  TEMP(K) =    57.45  PRESS =      0.00
 Etot   =       3.5310  EKtot   =      26.8870  EPtot      =     -23.3561
 BOND   =       6.2994  ANGLE   =      15.3763  DIHED      =      12.9274
 1-4 NB =       5.1100  1-4 EEL =      19.6567  VDWAALS    =      -2.3780
 EELEC  =     -80.3479  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


 NSTEP =    10  TIME(PS) =    0.010  TEMP(K) =    61.02  PRESS =      0.00
 Etot   =       8.4906  EKtot   =      28.5558  EPtot      =     -20.0652
 BOND   =       9.0485  ANGLE   =      18.3819  DIHED      =      11.6992
 1-4 NB =       5.3380  1-4 EEL =      16.4160  VDWAALS    =      -2.4726
 EELEC  =     -78.4762  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    20  TIME(PS) =    0.020  TEMP(K) =    71.90  PRESS =      0.00
 Etot   =      13.9648  EKtot   =      33.6469  EPtot      =     -19.6820
 BOND   =      13.2020  ANGLE   =      16.1126  DIHED      =       9.3872
 1-4 NB =       4.4568  1-4 EEL =      18.5683  VDWAALS    =      -0.7415
 EELEC  =     -80.6674  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    30  TIME(PS) =    0.030  TEMP(K) =    76.79  PRESS =      0.00
 Etot   =      19.1803  EKtot   =      35.9355  EPtot      =     -16.7552
 BOND   =      13.3291  ANGLE   =      17.0672  DIHED      =       8.9418
 1-4 NB =       4.7092  1-4 EEL =      23.5562  VDWAALS    =      -2.3706
 EELEC  =     -81.9883  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    40  TIME(PS) =    0.040  TEMP(K) =   101.28  PRESS =      0.00
 Etot   =      24.0438  EKtot   =      47.3996  EPtot      =     -23.3558
 BOND   =      10.2920  ANGLE   =      14.9634  DIHED      =       7.9428
 1-4 NB =       4.4907  1-4 EEL =      24.1672  VDWAALS    =      -2.1639
 EELEC  =     -83.0481  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    50  TIME(PS) =    0.050  TEMP(K) =    91.71  PRESS =      0.00
 Etot   =      28.8516  EKtot   =      42.9201  EPtot      =     -14.0685
 BOND   =       9.2663  ANGLE   =      23.2638  DIHED      =       8.0010
 1-4 NB =       6.8608  1-4 EEL =      23.1447  VDWAALS    =      -1.0064
 EELEC  =     -83.5988  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    60  TIME(PS) =    0.060  TEMP(K) =    90.40  PRESS =      0.00
 Etot   =      33.6243  EKtot   =      42.3046  EPtot      =      -8.6803
 BOND   =      12.2302  ANGLE   =      20.8793  DIHED      =       7.3831
 1-4 NB =       8.9735  1-4 EEL =      27.5591  VDWAALS    =      -0.8976
 EELEC  =     -84.8079  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    70  TIME(PS) =    0.070  TEMP(K) =   109.03  PRESS =      0.00
 Etot   =      38.2838  EKtot   =      51.0254  EPtot      =     -12.7416
 BOND   =      13.1860  ANGLE   =      17.8762  DIHED      =       9.2426
 1-4 NB =       7.2853  1-4 EEL =      27.4413  VDWAALS    =      -0.4866
 EELEC  =     -87.2865  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    80  TIME(PS) =    0.080  TEMP(K) =   113.14  PRESS =      0.00
 Etot   =      42.8582  EKtot   =      52.9477  EPtot      =     -10.0895
 BOND   =      13.1478  ANGLE   =      23.8640  DIHED      =       9.9536
 1-4 NB =       6.8597  1-4 EEL =      24.9833  VDWAALS    =      -0.4035
 EELEC  =     -88.4944  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =    90  TIME(PS) =    0.090  TEMP(K) =   104.79  PRESS =      0.00
 Etot   =      47.0872  EKtot   =      49.0382  EPtot      =      -1.9510
 BOND   =      17.7497  ANGLE   =      23.1455  DIHED      =       9.1027
 1-4 NB =       8.1694  1-4 EEL =      27.2788  VDWAALS    =      -0.1278
 EELEC  =     -87.2694  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   100  TIME(PS) =    0.100  TEMP(K) =   121.41  PRESS =      0.00
 Etot   =      51.2038  EKtot   =      56.8198  EPtot      =      -5.6160
 BOND   =      12.3809  ANGLE   =      23.6380  DIHED      =      10.1109
 1-4 NB =       6.3757  1-4 EEL =      26.4112  VDWAALS    =       0.0155
 EELEC  =     -84.5481  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.480399     ETOT(AT X=0) =  4.304E+00
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   110  TIME(PS) =    0.110  TEMP(K) =   126.35  PRESS =      0.00
 Etot   =      55.7117  EKtot   =      59.1311  EPtot      =      -3.4194
 BOND   =      19.2301  ANGLE   =      18.0740  DIHED      =      10.3957
 1-4 NB =       5.8505  1-4 EEL =      26.7594  VDWAALS    =      -1.0291
 EELEC  =     -82.7000  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   120  TIME(PS) =    0.120  TEMP(K) =   121.43  PRESS =      0.00
 Etot   =      59.7410  EKtot   =      56.8257  EPtot      =       2.9153
 BOND   =      14.1458  ANGLE   =      25.9106  DIHED      =      14.4835
 1-4 NB =       7.7421  1-4 EEL =      27.7902  VDWAALS    =      -1.7748
 EELEC  =     -85.3822  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   130  TIME(PS) =    0.130  TEMP(K) =   113.47  PRESS =      0.00
 Etot   =      63.5808  EKtot   =      53.1007  EPtot      =      10.4801
 BOND   =      11.4357  ANGLE   =      37.6059  DIHED      =      14.4597
 1-4 NB =       8.0316  1-4 EEL =      25.3450  VDWAALS    =      -1.8778
 EELEC  =     -84.5200  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   140  TIME(PS) =    0.140  TEMP(K) =   121.85  PRESS =      0.00
 Etot   =      67.6815  EKtot   =      57.0241  EPtot      =      10.6574
 BOND   =      11.6191  ANGLE   =      31.3995  DIHED      =      17.3472
 1-4 NB =       9.7084  1-4 EEL =      26.0457  VDWAALS    =      -2.6344
 EELEC  =     -82.8280  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   150  TIME(PS) =    0.150  TEMP(K) =   137.23  PRESS =      0.00
 Etot   =      71.8106  EKtot   =      64.2195  EPtot      =       7.5911
 BOND   =      15.4244  ANGLE   =      28.3535  DIHED      =      16.3149
 1-4 NB =       5.6734  1-4 EEL =      26.4789  VDWAALS    =      -2.6104
 EELEC  =     -82.0436  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   160  TIME(PS) =    0.160  TEMP(K) =   132.14  PRESS =      0.00
 Etot   =      75.9306  EKtot   =      61.8410  EPtot      =      14.0896
 BOND   =      15.1082  ANGLE   =      36.6210  DIHED      =      19.3551
 1-4 NB =       5.4543  1-4 EEL =      23.0216  VDWAALS    =      -3.6888
 EELEC  =     -81.7818  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   170  TIME(PS) =    0.170  TEMP(K) =   109.42  PRESS =      0.00
 Etot   =      79.7779  EKtot   =      51.2080  EPtot      =      28.5699
 BOND   =      17.0698  ANGLE   =      42.0966  DIHED      =      22.1244
 1-4 NB =       7.9254  1-4 EEL =      25.8846  VDWAALS    =      -1.0759
 EELEC  =     -85.4550  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   180  TIME(PS) =    0.180  TEMP(K) =   129.77  PRESS =      0.00
 Etot   =      83.1417  EKtot   =      60.7322  EPtot      =      22.4095
 BOND   =      18.0007  ANGLE   =      40.5457  DIHED      =      20.5305
 1-4 NB =       4.3053  1-4 EEL =      32.1589  VDWAALS    =      -3.1934
 EELEC  =     -89.9382  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   190  TIME(PS) =    0.190  TEMP(K) =   142.87  PRESS =      0.00
 Etot   =      86.8709  EKtot   =      66.8613  EPtot      =      20.0096
 BOND   =      25.8831  ANGLE   =      30.0472  DIHED      =      15.5235
 1-4 NB =       5.3594  1-4 EEL =      31.3811  VDWAALS    =      -1.1382
 EELEC  =     -87.0466  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   200  TIME(PS) =    0.200  TEMP(K) =   156.69  PRESS =      0.00
 Etot   =      90.5329  EKtot   =      73.3267  EPtot      =      17.2062
 BOND   =      18.0883  ANGLE   =      30.3655  DIHED      =      22.8623
 1-4 NB =       6.6875  1-4 EEL =      27.0722  VDWAALS    =      -2.3057
 EELEC  =     -85.5639  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.395712     ETOT(AT X=0) =  5.169E+01
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   210  TIME(PS) =    0.210  TEMP(K) =   176.35  PRESS =      0.00
 Etot   =      94.1226  EKtot   =      82.5306  EPtot      =      11.5920
 BOND   =      18.5993  ANGLE   =      25.7664  DIHED      =      18.0745
 1-4 NB =       9.5029  1-4 EEL =      27.8464  VDWAALS    =      -1.4890
 EELEC  =     -86.7085  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   220  TIME(PS) =    0.220  TEMP(K) =   172.67  PRESS =      0.00
 Etot   =      96.8993  EKtot   =      80.8088  EPtot      =      16.0904
 BOND   =      18.3704  ANGLE   =      33.3965  DIHED      =      14.7391
 1-4 NB =      11.7004  1-4 EEL =      29.8116  VDWAALS    =      -2.9244
 EELEC  =     -89.0032  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   230  TIME(PS) =    0.230  TEMP(K) =   187.41  PRESS =      0.00
 Etot   =      99.3179  EKtot   =      87.7069  EPtot      =      11.6109
 BOND   =      15.1000  ANGLE   =      31.0027  DIHED      =      15.4290
 1-4 NB =      10.5388  1-4 EEL =      31.7177  VDWAALS    =      -0.9042
 EELEC  =     -91.2731  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   240  TIME(PS) =    0.240  TEMP(K) =   187.66  PRESS =      0.00
 Etot   =     102.1289  EKtot   =      87.8211  EPtot      =      14.3078
 BOND   =      19.0458  ANGLE   =      34.1657  DIHED      =      13.6275
 1-4 NB =      11.5560  1-4 EEL =      27.6867  VDWAALS    =      -2.3563
 EELEC  =     -89.4174  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   250  TIME(PS) =    0.250  TEMP(K) =   173.04  PRESS =      0.00
 Etot   =     105.4501  EKtot   =      80.9817  EPtot      =      24.4684
 BOND   =      24.3035  ANGLE   =      37.7658  DIHED      =      18.7652
 1-4 NB =       8.3693  1-4 EEL =      24.8493  VDWAALS    =      -1.0940
 EELEC  =     -88.4906  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   260  TIME(PS) =    0.260  TEMP(K) =   208.46  PRESS =      0.00
 Etot   =     107.4072  EKtot   =      97.5579  EPtot      =       9.8493
 BOND   =      15.1644  ANGLE   =      32.7645  DIHED      =      17.7395
 1-4 NB =       6.8858  1-4 EEL =      28.3260  VDWAALS    =      -2.1842
 EELEC  =     -88.8466  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   270  TIME(PS) =    0.270  TEMP(K) =   212.55  PRESS =      0.00
 Etot   =     109.8696  EKtot   =      99.4700  EPtot      =      10.3996
 BOND   =      18.4559  ANGLE   =      28.5431  DIHED      =      20.7129
 1-4 NB =       4.5561  1-4 EEL =      29.4268  VDWAALS    =      -1.1199
 EELEC  =     -90.1752  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   280  TIME(PS) =    0.280  TEMP(K) =   169.62  PRESS =      0.00
 Etot   =     112.9387  EKtot   =      79.3786  EPtot      =      33.5601
 BOND   =      34.8307  ANGLE   =      31.4250  DIHED      =      19.8493
 1-4 NB =       5.5075  1-4 EEL =      26.5507  VDWAALS    =       2.6440
 EELEC  =     -87.2472  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   290  TIME(PS) =    0.290  TEMP(K) =   168.75  PRESS =      0.00
 Etot   =     115.6922  EKtot   =      78.9712  EPtot      =      36.7210
 BOND   =      26.2396  ANGLE   =      40.4867  DIHED      =      20.8652
 1-4 NB =       6.1871  1-4 EEL =      24.2993  VDWAALS    =       3.5513
 EELEC  =     -84.9081  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   300  TIME(PS) =    0.300  TEMP(K) =   177.92  PRESS =      0.00
 Etot   =     118.7216  EKtot   =      83.2634  EPtot      =      35.4582
 BOND   =      26.3677  ANGLE   =      38.0525  DIHED      =      19.1205
 1-4 NB =       5.9646  1-4 EEL =      25.5370  VDWAALS    =       4.5324
 EELEC  =     -84.1165  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.271657     ETOT(AT X=0) =  9.137E+01
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   310  TIME(PS) =    0.310  TEMP(K) =   189.74  PRESS =      0.00
 Etot   =     121.4939  EKtot   =      88.7953  EPtot      =      32.6986
 BOND   =      28.9944  ANGLE   =      39.2157  DIHED      =      16.7254
 1-4 NB =       6.3006  1-4 EEL =      22.9555  VDWAALS    =       1.3758
 EELEC  =     -82.8688  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   320  TIME(PS) =    0.320  TEMP(K) =   186.93  PRESS =      0.00
 Etot   =     123.4620  EKtot   =      87.4791  EPtot      =      35.9829
 BOND   =      25.0301  ANGLE   =      48.7578  DIHED      =      19.2168
 1-4 NB =       7.8789  1-4 EEL =      18.0210  VDWAALS    =      -0.4631
 EELEC  =     -82.4586  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   330  TIME(PS) =    0.330  TEMP(K) =   224.89  PRESS =      0.00
 Etot   =     126.1405  EKtot   =     105.2442  EPtot      =      20.8963
 BOND   =      25.9173  ANGLE   =      38.7418  DIHED      =      13.4751
 1-4 NB =       5.4409  1-4 EEL =      21.4125  VDWAALS    =      -1.9337
 EELEC  =     -82.1576  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   340  TIME(PS) =    0.340  TEMP(K) =   229.10  PRESS =      0.00
 Etot   =     129.1833  EKtot   =     107.2156  EPtot      =      21.9677
 BOND   =      28.7214  ANGLE   =      36.5901  DIHED      =      13.0123
 1-4 NB =       4.5133  1-4 EEL =      25.4720  VDWAALS    =      -2.4254
 EELEC  =     -83.9160  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   350  TIME(PS) =    0.350  TEMP(K) =   202.42  PRESS =      0.00
 Etot   =     131.4644  EKtot   =      94.7309  EPtot      =      36.7335
 BOND   =      31.5150  ANGLE   =      49.2389  DIHED      =      13.7777
 1-4 NB =       5.8825  1-4 EEL =      24.1695  VDWAALS    =      -1.5456
 EELEC  =     -86.3045  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   360  TIME(PS) =    0.360  TEMP(K) =   187.72  PRESS =      0.00
 Etot   =     132.8221  EKtot   =      87.8480  EPtot      =      44.9741
 BOND   =      27.8419  ANGLE   =      51.1419  DIHED      =      14.4308
 1-4 NB =       9.9999  1-4 EEL =      24.8920  VDWAALS    =       3.2594
 EELEC  =     -86.5919  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   370  TIME(PS) =    0.370  TEMP(K) =   195.11  PRESS =      0.00
 Etot   =     134.2850  EKtot   =      91.3105  EPtot      =      42.9745
 BOND   =      20.3323  ANGLE   =      58.3425  DIHED      =      17.0215
 1-4 NB =      10.2622  1-4 EEL =      25.3581  VDWAALS    =       0.5707
 EELEC  =     -88.9127  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   380  TIME(PS) =    0.380  TEMP(K) =   243.82  PRESS =      0.00
 Etot   =     135.4115  EKtot   =     114.1036  EPtot      =      21.3079
 BOND   =      23.7529  ANGLE   =      35.5289  DIHED      =      12.4747
 1-4 NB =       8.3672  1-4 EEL =      23.1226  VDWAALS    =       0.2120
 EELEC  =     -82.1504  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   390  TIME(PS) =    0.390  TEMP(K) =   231.24  PRESS =      0.00
 Etot   =     137.5653  EKtot   =     108.2175  EPtot      =      29.3478
 BOND   =      25.6388  ANGLE   =      39.0795  DIHED      =      14.8487
 1-4 NB =      10.7551  1-4 EEL =      20.0970  VDWAALS    =      -1.2021
 EELEC  =     -79.8691  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   400  TIME(PS) =    0.400  TEMP(K) =   201.91  PRESS =      0.00
 Etot   =     140.5702  EKtot   =      94.4909  EPtot      =      46.0792
 BOND   =      36.4196  ANGLE   =      40.9665  DIHED      =      18.8709
 1-4 NB =      10.6672  1-4 EEL =      14.9307  VDWAALS    =       1.6484
 EELEC  =     -77.4241  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.204642     ETOT(AT X=0) =  1.200E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   410  TIME(PS) =    0.410  TEMP(K) =   197.90  PRESS =      0.00
 Etot   =     142.9905  EKtot   =      92.6131  EPtot      =      50.3773
 BOND   =      34.3420  ANGLE   =      52.3830  DIHED      =      17.2314
 1-4 NB =       6.2616  1-4 EEL =      20.8379  VDWAALS    =      -0.1256
 EELEC  =     -80.5529  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   420  TIME(PS) =    0.420  TEMP(K) =   231.17  PRESS =      0.00
 Etot   =     143.8212  EKtot   =     108.1858  EPtot      =      35.6354
 BOND   =      32.1175  ANGLE   =      38.3104  DIHED      =      15.8833
 1-4 NB =       6.0869  1-4 EEL =      23.6144  VDWAALS    =      -1.8077
 EELEC  =     -78.5695  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   430  TIME(PS) =    0.430  TEMP(K) =   231.53  PRESS =      0.00
 Etot   =     145.2992  EKtot   =     108.3520  EPtot      =      36.9472
 BOND   =      32.9901  ANGLE   =      36.1720  DIHED      =      20.6963
 1-4 NB =       4.4887  1-4 EEL =      23.1612  VDWAALS    =      -0.7618
 EELEC  =     -79.7993  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   440  TIME(PS) =    0.440  TEMP(K) =   181.66  PRESS =      0.00
 Etot   =     147.7164  EKtot   =      85.0147  EPtot      =      62.7016
 BOND   =      36.2524  ANGLE   =      57.3760  DIHED      =      22.2476
 1-4 NB =       6.2869  1-4 EEL =      17.8559  VDWAALS    =      -2.2436
 EELEC  =     -75.0736  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   450  TIME(PS) =    0.450  TEMP(K) =   199.37  PRESS =      0.00
 Etot   =     149.4753  EKtot   =      93.2999  EPtot      =      56.1753
 BOND   =      41.2545  ANGLE   =      40.8744  DIHED      =      17.4355
 1-4 NB =       7.5088  1-4 EEL =      23.1774  VDWAALS    =       3.8772
 EELEC  =     -77.9524  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   460  TIME(PS) =    0.460  TEMP(K) =   211.22  PRESS =      0.00
 Etot   =     151.0491  EKtot   =      98.8484  EPtot      =      52.2007
 BOND   =      39.2373  ANGLE   =      40.7723  DIHED      =      22.6235
 1-4 NB =       7.5210  1-4 EEL =      22.4158  VDWAALS    =      -2.0485
 EELEC  =     -78.3208  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   470  TIME(PS) =    0.470  TEMP(K) =   242.61  PRESS =      0.00
 Etot   =     152.9246  EKtot   =     113.5364  EPtot      =      39.3882
 BOND   =      26.6351  ANGLE   =      45.2889  DIHED      =      18.2006
 1-4 NB =       6.7236  1-4 EEL =      16.9788  VDWAALS    =      -0.4578
 EELEC  =     -73.9809  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   480  TIME(PS) =    0.480  TEMP(K) =   209.77  PRESS =      0.00
 Etot   =     154.9767  EKtot   =      98.1668  EPtot      =      56.8099
 BOND   =      43.0649  ANGLE   =      50.2243  DIHED      =      16.8303
 1-4 NB =       6.6687  1-4 EEL =      13.4674  VDWAALS    =      -2.5518
 EELEC  =     -70.8940  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   490  TIME(PS) =    0.490  TEMP(K) =   233.10  PRESS =      0.00
 Etot   =     156.3771  EKtot   =     109.0879  EPtot      =      47.2892
 BOND   =      31.2255  ANGLE   =      50.6639  DIHED      =      19.9999
 1-4 NB =       7.4482  1-4 EEL =      11.2483  VDWAALS    =      -2.7592
 EELEC  =     -70.5373  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   500  TIME(PS) =    0.500  TEMP(K) =   252.22  PRESS =      0.00
 Etot   =     157.0948  EKtot   =     118.0358  EPtot      =      39.0590
 BOND   =      32.7951  ANGLE   =      45.1362  DIHED      =      12.6139
 1-4 NB =       6.6400  1-4 EEL =      18.5031  VDWAALS    =      -2.3583
 EELEC  =     -74.2710  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.175983     ETOT(AT X=0) =  1.403E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   510  TIME(PS) =    0.510  TEMP(K) =   272.42  PRESS =      0.00
 Etot   =     157.7360  EKtot   =     127.4898  EPtot      =      30.2462
 BOND   =      20.8557  ANGLE   =      42.8819  DIHED      =      11.6733
 1-4 NB =      11.3181  1-4 EEL =      21.0905  VDWAALS    =       1.2715
 EELEC  =     -78.8449  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   520  TIME(PS) =    0.520  TEMP(K) =   239.85  PRESS =      0.00
 Etot   =     159.7830  EKtot   =     112.2466  EPtot      =      47.5364
 BOND   =      42.8500  ANGLE   =      39.7372  DIHED      =      22.6024
 1-4 NB =       6.2520  1-4 EEL =      21.1010  VDWAALS    =      -0.6189
 EELEC  =     -84.3873  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   530  TIME(PS) =    0.530  TEMP(K) =   253.54  PRESS =      0.00
 Etot   =     161.5409  EKtot   =     118.6530  EPtot      =      42.8879
 BOND   =      35.3413  ANGLE   =      37.5055  DIHED      =      21.8272
 1-4 NB =       8.2647  1-4 EEL =      29.2086  VDWAALS    =      -0.9454
 EELEC  =     -88.3139  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   540  TIME(PS) =    0.540  TEMP(K) =   228.67  PRESS =      0.00
 Etot   =     161.9053  EKtot   =     107.0139  EPtot      =      54.8914
 BOND   =      32.9290  ANGLE   =      53.0285  DIHED      =      13.2752
 1-4 NB =      10.7123  1-4 EEL =      38.3671  VDWAALS    =      -1.5830
 EELEC  =     -91.8377  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   550  TIME(PS) =    0.550  TEMP(K) =   264.07  PRESS =      0.00
 Etot   =     162.0536  EKtot   =     123.5815  EPtot      =      38.4722
 BOND   =      27.3307  ANGLE   =      43.1483  DIHED      =      18.3735
 1-4 NB =       9.2033  1-4 EEL =      32.8580  VDWAALS    =      -0.4330
 EELEC  =     -92.0086  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   560  TIME(PS) =    0.560  TEMP(K) =   257.03  PRESS =      0.00
 Etot   =     163.1763  EKtot   =     120.2856  EPtot      =      42.8908
 BOND   =      24.0202  ANGLE   =      54.8902  DIHED      =      18.2798
 1-4 NB =       6.8973  1-4 EEL =      29.0589  VDWAALS    =      -1.7917
 EELEC  =     -88.4639  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   570  TIME(PS) =    0.570  TEMP(K) =   275.33  PRESS =      0.00
 Etot   =     164.3026  EKtot   =     128.8491  EPtot      =      35.4535
 BOND   =      24.1194  ANGLE   =      42.4999  DIHED      =      17.5589
 1-4 NB =       9.4380  1-4 EEL =      29.1234  VDWAALS    =      -1.6939
 EELEC  =     -85.5922  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   580  TIME(PS) =    0.580  TEMP(K) =   207.54  PRESS =      0.00
 Etot   =     166.7960  EKtot   =      97.1240  EPtot      =      69.6720
 BOND   =      34.8496  ANGLE   =      57.6338  DIHED      =      29.2061
 1-4 NB =       8.9078  1-4 EEL =      22.0881  VDWAALS    =      -0.5265
 EELEC  =     -82.4869  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   590  TIME(PS) =    0.590  TEMP(K) =   211.51  PRESS =      0.00
 Etot   =     168.4149  EKtot   =      98.9854  EPtot      =      69.4295
 BOND   =      28.3042  ANGLE   =      59.9120  DIHED      =      36.0079
 1-4 NB =       9.4448  1-4 EEL =      19.0353  VDWAALS    =      -1.3448
 EELEC  =     -81.9298  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   600  TIME(PS) =    0.600  TEMP(K) =   251.83  PRESS =      0.00
 Etot   =     169.9990  EKtot   =     117.8533  EPtot      =      52.1457
 BOND   =      39.7815  ANGLE   =      33.3190  DIHED      =      30.5192
 1-4 NB =       8.2189  1-4 EEL =      21.1436  VDWAALS    =      -0.7231
 EELEC  =     -80.1135  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.119837     ETOT(AT X=0) =  1.570E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   610  TIME(PS) =    0.610  TEMP(K) =   278.36  PRESS =      0.00
 Etot   =     170.9755  EKtot   =     130.2661  EPtot      =      40.7094
 BOND   =      28.2258  ANGLE   =      40.8697  DIHED      =      24.6941
 1-4 NB =       5.1351  1-4 EEL =      25.0392  VDWAALS    =      -1.7361
 EELEC  =     -81.5184  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   620  TIME(PS) =    0.620  TEMP(K) =   275.48  PRESS =      0.00
 Etot   =     171.5657  EKtot   =     128.9205  EPtot      =      42.6453
 BOND   =      28.8897  ANGLE   =      42.7165  DIHED      =      25.4957
 1-4 NB =       6.7026  1-4 EEL =      24.2184  VDWAALS    =      -2.5114
 EELEC  =     -82.8663  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   630  TIME(PS) =    0.630  TEMP(K) =   255.78  PRESS =      0.00
 Etot   =     172.8146  EKtot   =     119.7021  EPtot      =      53.1124
 BOND   =      33.7339  ANGLE   =      47.7918  DIHED      =      20.7514
 1-4 NB =      10.6312  1-4 EEL =      25.1014  VDWAALS    =       0.0481
 EELEC  =     -84.9454  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   640  TIME(PS) =    0.640  TEMP(K) =   263.89  PRESS =      0.00
 Etot   =     173.1232  EKtot   =     123.4964  EPtot      =      49.6268
 BOND   =      27.0266  ANGLE   =      52.8342  DIHED      =      21.5278
 1-4 NB =       9.1250  1-4 EEL =      24.2856  VDWAALS    =      -1.3979
 EELEC  =     -83.7743  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   650  TIME(PS) =    0.650  TEMP(K) =   274.59  PRESS =      0.00
 Etot   =     173.5012  EKtot   =     128.5048  EPtot      =      44.9964
 BOND   =      22.7203  ANGLE   =      39.2780  DIHED      =      31.4339
 1-4 NB =      10.0433  1-4 EEL =      22.0879  VDWAALS    =      -1.3852
 EELEC  =     -79.1819  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   660  TIME(PS) =    0.660  TEMP(K) =   244.44  PRESS =      0.00
 Etot   =     176.1097  EKtot   =     114.3917  EPtot      =      61.7180
 BOND   =      40.4384  ANGLE   =      46.3789  DIHED      =      26.6254
 1-4 NB =       8.0833  1-4 EEL =      20.0325  VDWAALS    =      -2.5217
 EELEC  =     -77.3188  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   670  TIME(PS) =    0.670  TEMP(K) =   315.83  PRESS =      0.00
 Etot   =     175.8483  EKtot   =     147.8037  EPtot      =      28.0446
 BOND   =      37.4492  ANGLE   =      29.4081  DIHED      =      11.2214
 1-4 NB =       7.1376  1-4 EEL =      21.0006  VDWAALS    =       0.1815
 EELEC  =     -78.3538  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   680  TIME(PS) =    0.680  TEMP(K) =   299.95  PRESS =      0.00
 Etot   =     175.3229  EKtot   =     140.3710  EPtot      =      34.9519
 BOND   =      32.5835  ANGLE   =      37.2182  DIHED      =      19.2016
 1-4 NB =       5.6676  1-4 EEL =      19.9106  VDWAALS    =      -2.5233
 EELEC  =     -77.1063  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   690  TIME(PS) =    0.690  TEMP(K) =   277.72  PRESS =      0.00
 Etot   =     176.4341  EKtot   =     129.9705  EPtot      =      46.4636
 BOND   =      35.6674  ANGLE   =      42.2111  DIHED      =      23.8084
 1-4 NB =       4.6895  1-4 EEL =      18.0331  VDWAALS    =      -2.2069
 EELEC  =     -75.7389  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   700  TIME(PS) =    0.700  TEMP(K) =   241.63  PRESS =      0.00
 Etot   =     177.8872  EKtot   =     113.0811  EPtot      =      64.8061
 BOND   =      36.6941  ANGLE   =      57.3182  DIHED      =      20.6611
 1-4 NB =       7.0350  1-4 EEL =      18.9189  VDWAALS    =      -1.0470
 EELEC  =     -74.7741  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.071714     ETOT(AT X=0) =  1.704E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   710  TIME(PS) =    0.710  TEMP(K) =   267.76  PRESS =      0.00
 Etot   =     178.2446  EKtot   =     125.3071  EPtot      =      52.9375
 BOND   =      31.7638  ANGLE   =      50.5973  DIHED      =      19.9497
 1-4 NB =       7.8166  1-4 EEL =      15.4619  VDWAALS    =      -0.1789
 EELEC  =     -72.4728  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   720  TIME(PS) =    0.720  TEMP(K) =   262.27  PRESS =      0.00
 Etot   =     178.5225  EKtot   =     122.7374  EPtot      =      55.7851
 BOND   =      46.0710  ANGLE   =      41.9927  DIHED      =      19.8356
 1-4 NB =       7.6089  1-4 EEL =      16.1886  VDWAALS    =      -2.3810
 EELEC  =     -73.5308  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   730  TIME(PS) =    0.730  TEMP(K) =   296.56  PRESS =      0.00
 Etot   =     177.7112  EKtot   =     138.7832  EPtot      =      38.9281
 BOND   =      27.1864  ANGLE   =      41.8337  DIHED      =      21.0011
 1-4 NB =       7.2979  1-4 EEL =      21.1715  VDWAALS    =      -1.9835
 EELEC  =     -77.5791  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   740  TIME(PS) =    0.740  TEMP(K) =   271.73  PRESS =      0.00
 Etot   =     178.8713  EKtot   =     127.1646  EPtot      =      51.7067
 BOND   =      30.4176  ANGLE   =      42.8800  DIHED      =      27.3852
 1-4 NB =       7.9030  1-4 EEL =      20.5692  VDWAALS    =      -1.0805
 EELEC  =     -76.3678  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   750  TIME(PS) =    0.750  TEMP(K) =   227.29  PRESS =      0.00
 Etot   =     180.9072  EKtot   =     106.3677  EPtot      =      74.5396
 BOND   =      38.4802  ANGLE   =      59.5697  DIHED      =      23.1857
 1-4 NB =      10.6426  1-4 EEL =      17.9211  VDWAALS    =      -0.3781
 EELEC  =     -74.8815  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   760  TIME(PS) =    0.760  TEMP(K) =   267.80  PRESS =      0.00
 Etot   =     181.5024  EKtot   =     125.3280  EPtot      =      56.1744
 BOND   =      31.0731  ANGLE   =      54.7456  DIHED      =      18.2842
 1-4 NB =       9.6250  1-4 EEL =      17.3641  VDWAALS    =      -1.1146
 EELEC  =     -73.8030  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   770  TIME(PS) =    0.770  TEMP(K) =   288.34  PRESS =      0.00
 Etot   =     181.6535  EKtot   =     134.9384  EPtot      =      46.7152
 BOND   =      27.9948  ANGLE   =      38.0618  DIHED      =      19.6436
 1-4 NB =       9.7406  1-4 EEL =      17.9557  VDWAALS    =       7.8407
 EELEC  =     -74.5220  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   780  TIME(PS) =    0.780  TEMP(K) =   250.55  PRESS =      0.00
 Etot   =     182.7572  EKtot   =     117.2527  EPtot      =      65.5045
 BOND   =      30.6756  ANGLE   =      61.6721  DIHED      =      20.1926
 1-4 NB =       9.2452  1-4 EEL =      14.1677  VDWAALS    =       1.1685
 EELEC  =     -71.6172  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   790  TIME(PS) =    0.790  TEMP(K) =   260.06  PRESS =      0.00
 Etot   =     184.1655  EKtot   =     121.7060  EPtot      =      62.4595
 BOND   =      36.2277  ANGLE   =      50.8653  DIHED      =      15.3703
 1-4 NB =      13.4078  1-4 EEL =      13.7163  VDWAALS    =       3.6593
 EELEC  =     -70.7872  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   800  TIME(PS) =    0.800  TEMP(K) =   272.40  PRESS =      0.00
 Etot   =     185.2770  EKtot   =     127.4776  EPtot      =      57.7995
 BOND   =      38.8336  ANGLE   =      52.0993  DIHED      =      17.4897
 1-4 NB =       8.9358  1-4 EEL =      16.4344  VDWAALS    =      -0.1249
 EELEC  =     -75.8684  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.082027     ETOT(AT X=0) =  1.765E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   810  TIME(PS) =    0.810  TEMP(K) =   278.30  PRESS =      0.00
 Etot   =     184.6091  EKtot   =     130.2415  EPtot      =      54.3676
 BOND   =      27.0605  ANGLE   =      60.9392  DIHED      =      18.4724
 1-4 NB =       3.7153  1-4 EEL =      22.4779  VDWAALS    =       1.6549
 EELEC  =     -79.9524  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   820  TIME(PS) =    0.820  TEMP(K) =   301.70  PRESS =      0.00
 Etot   =     185.1561  EKtot   =     141.1908  EPtot      =      43.9652
 BOND   =      31.9360  ANGLE   =      46.5233  DIHED      =      16.6140
 1-4 NB =       4.2610  1-4 EEL =      20.2255  VDWAALS    =       1.5096
 EELEC  =     -77.1041  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   830  TIME(PS) =    0.830  TEMP(K) =   258.90  PRESS =      0.00
 Etot   =     186.6502  EKtot   =     121.1617  EPtot      =      65.4886
 BOND   =      45.9803  ANGLE   =      61.8185  DIHED      =      11.7896
 1-4 NB =       4.7686  1-4 EEL =      19.2863  VDWAALS    =       0.3975
 EELEC  =     -78.5523  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   840  TIME(PS) =    0.840  TEMP(K) =   298.93  PRESS =      0.00
 Etot   =     187.9144  EKtot   =     139.8965  EPtot      =      48.0179
 BOND   =      46.9488  ANGLE   =      39.2684  DIHED      =       9.5938
 1-4 NB =      10.7107  1-4 EEL =      19.0507  VDWAALS    =      -1.4830
 EELEC  =     -76.0715  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   850  TIME(PS) =    0.850  TEMP(K) =   295.59  PRESS =      0.00
 Etot   =     187.6527  EKtot   =     138.3319  EPtot      =      49.3208
 BOND   =      35.1650  ANGLE   =      53.1692  DIHED      =       9.3635
 1-4 NB =      12.0939  1-4 EEL =      18.8099  VDWAALS    =      -2.0290
 EELEC  =     -77.2516  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   860  TIME(PS) =    0.860  TEMP(K) =   306.44  PRESS =      0.00
 Etot   =     186.6760  EKtot   =     143.4093  EPtot      =      43.2667
 BOND   =      37.1619  ANGLE   =      43.2511  DIHED      =      12.4826
 1-4 NB =      11.5349  1-4 EEL =      16.6208  VDWAALS    =       0.3452
 EELEC  =     -78.1298  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   870  TIME(PS) =    0.870  TEMP(K) =   300.04  PRESS =      0.00
 Etot   =     186.9415  EKtot   =     140.4156  EPtot      =      46.5258
 BOND   =      45.2209  ANGLE   =      36.8158  DIHED      =      13.2959
 1-4 NB =      10.1246  1-4 EEL =      17.6588  VDWAALS    =       0.9207
 EELEC  =     -77.5107  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   880  TIME(PS) =    0.880  TEMP(K) =   257.00  PRESS =      0.00
 Etot   =     187.5455  EKtot   =     120.2728  EPtot      =      67.2727
 BOND   =      49.2860  ANGLE   =      56.9846  DIHED      =      10.4027
 1-4 NB =       9.4618  1-4 EEL =      14.6212  VDWAALS    =       1.2970
 EELEC  =     -74.7806  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   890  TIME(PS) =    0.890  TEMP(K) =   275.08  PRESS =      0.00
 Etot   =     189.2708  EKtot   =     128.7340  EPtot      =      60.5368
 BOND   =      49.8972  ANGLE   =      51.1391  DIHED      =      12.8982
 1-4 NB =       6.5422  1-4 EEL =      14.1148  VDWAALS    =      -1.0626
 EELEC  =     -72.9921  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   900  TIME(PS) =    0.900  TEMP(K) =   268.67  PRESS =      0.00
 Etot   =     189.3012  EKtot   =     125.7340  EPtot      =      63.5672
 BOND   =      51.6183  ANGLE   =      50.8299  DIHED      =      14.3294
 1-4 NB =       5.1381  1-4 EEL =      17.2598  VDWAALS    =      -0.1315
 EELEC  =     -75.4769  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.033178     ETOT(AT X=0) =  1.854E+02
 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   910  TIME(PS) =    0.910  TEMP(K) =   300.69  PRESS =      0.00
 Etot   =     189.2802  EKtot   =     140.7188  EPtot      =      48.5614
 BOND   =      46.9254  ANGLE   =      44.9707  DIHED      =      16.0654
 1-4 NB =       4.9749  1-4 EEL =      15.0313  VDWAALS    =       0.3402
 EELEC  =     -79.7465  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   920  TIME(PS) =    0.920  TEMP(K) =   298.04  PRESS =      0.00
 Etot   =     189.5456  EKtot   =     139.4758  EPtot      =      50.0698
 BOND   =      32.2665  ANGLE   =      55.4807  DIHED      =      13.0095
 1-4 NB =       9.3936  1-4 EEL =      18.6664  VDWAALS    =       4.8782
 EELEC  =     -83.6249  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   930  TIME(PS) =    0.930  TEMP(K) =   249.40  PRESS =      0.00
 Etot   =     190.3035  EKtot   =     116.7137  EPtot      =      73.5898
 BOND   =      39.1451  ANGLE   =      64.3348  DIHED      =      16.9290
 1-4 NB =       9.6350  1-4 EEL =      26.7247  VDWAALS    =       3.1267
 EELEC  =     -86.3056  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   940  TIME(PS) =    0.940  TEMP(K) =   278.13  PRESS =      0.00
 Etot   =     190.5001  EKtot   =     130.1621  EPtot      =      60.3380
 BOND   =      49.0087  ANGLE   =      41.1371  DIHED      =      20.5595
 1-4 NB =       7.9398  1-4 EEL =      26.5729  VDWAALS    =       7.2998
 EELEC  =     -92.1798  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   950  TIME(PS) =    0.950  TEMP(K) =   242.43  PRESS =      0.00
 Etot   =     192.3355  EKtot   =     113.4513  EPtot      =      78.8842
 BOND   =      39.5771  ANGLE   =      64.4123  DIHED      =      23.7025
 1-4 NB =      12.0724  1-4 EEL =      23.8817  VDWAALS    =       4.9091
 EELEC  =     -89.6710  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   960  TIME(PS) =    0.960  TEMP(K) =   297.32  PRESS =      0.00
 Etot   =     192.1173  EKtot   =     139.1388  EPtot      =      52.9785
 BOND   =      34.9521  ANGLE   =      49.5058  DIHED      =      20.9578
 1-4 NB =       6.1798  1-4 EEL =      28.9694  VDWAALS    =       0.9582
 EELEC  =     -88.5446  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   970  TIME(PS) =    0.970  TEMP(K) =   304.95  PRESS =      0.00
 Etot   =     193.2883  EKtot   =     142.7108  EPtot      =      50.5775
 BOND   =      45.7788  ANGLE   =      43.9144  DIHED      =      11.8462
 1-4 NB =       5.5216  1-4 EEL =      22.8974  VDWAALS    =       2.0165
 EELEC  =     -81.3975  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   980  TIME(PS) =    0.980  TEMP(K) =   256.83  PRESS =      0.00
 Etot   =     194.1619  EKtot   =     120.1900  EPtot      =      73.9719
 BOND   =      47.4913  ANGLE   =      60.1195  DIHED      =      21.7069
 1-4 NB =       5.4460  1-4 EEL =      14.7750  VDWAALS    =      -0.3331
 EELEC  =     -75.2337  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =   990  TIME(PS) =    0.990  TEMP(K) =   268.63  PRESS =      0.00
 Etot   =     194.2782  EKtot   =     125.7148  EPtot      =      68.5634
 BOND   =      56.5911  ANGLE   =      47.1953  DIHED      =      20.1009
 1-4 NB =       5.7464  1-4 EEL =      10.1853  VDWAALS    =       1.8201
 EELEC  =     -73.0757  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

 NB-update: NPAIRS =    2391  HBPAIR =       0

 NSTEP =  1000  TIME(PS) =    1.000  TEMP(K) =   297.06  PRESS =      0.00
 Etot   =     192.7903  EKtot   =     139.0203  EPtot      =      53.7700
 BOND   =      22.4020  ANGLE   =      66.6662  DIHED      =      18.3047
 1-4 NB =       7.5314  1-4 EEL =      10.4345  VDWAALS    =      -1.4530
 EELEC  =     -70.1158  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      RESULT OF LEAST SQUARE FIT OVER  100 STEPS
      ENERGY DRIFT PER STEP =      0.055938     ETOT(AT X=0) =  1.888E+02

      A V E R A G E S   O V E R  1000 S T E P S


 NSTEP =  1000  TIME(PS) =    1.000  TEMP(K) =   216.75  PRESS =      0.00
 Etot   =     138.1320  EKtot   =     101.4357  EPtot      =      36.6963
 BOND   =      29.2846  ANGLE   =      41.8705  DIHED      =      17.3728
 1-4 NB =       7.6251  1-4 EEL =      22.3503  VDWAALS    =      -0.3453
 EELEC  =     -81.4616  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------


      R M S  F L U C T U A T I O N S


 NSTEP =  1000  TIME(PS) =    1.000  TEMP(K) =    65.05  PRESS =      0.00
 Etot   =      52.1032  EKtot   =      30.4416  EPtot      =      24.8429
 BOND   =      10.9344  ANGLE   =      12.3544  DIHED      =       5.3952
 1-4 NB =       2.2233  1-4 EEL =       5.2386  VDWAALS    =       2.2405
 EELEC  =       5.5715  EHBOND  =       0.0000  CONSTRAINT =       0.0000
 ------------------------------------------------------------------------------

|         ELAPSED TIME =     13.390     TOTAL TIME =     13.390

      Routine         Sec       %
      ----------------------------
|     Pairlist       0.21    1.57
|     Nonbond        1.65   12.32
|     Bond           0.14    1.05
|     Angle          1.55   11.58
|     Dihedral       7.86   58.70
|     Shake          0.00    0.00
|     Quick3         0.00    0.00
|     Force          0.00    0.00
|     Other          1.98   14.79
      ----------------------------
|     Total         13.39    0.00 Hours

|     Nonsetup      13.14   98.13%

|     Setup wallclock           1 seconds
|     Nonsetup wallclock       17 seconds
