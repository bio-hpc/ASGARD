C
C     PM3 SEMIEMPIRICAL PARAMETERS:
C     ----------------------------
C
C     DATA FOR HYDROGEN:
C
      DATA UCORE3(0,1)    /  -13.0733210D0/
      DATA EXPNT3(0,1)    /    0.9678070D0/
      DATA AL3(0,1)       /    0.5437048D0/

C     note that the following 2 entries didn't exist
C     in Steve's original block.F:
      DATA AL3(1,1)       /    0.5437048D0/
      DATA AL3(2,1)       /    0.5437048D0/

      DATA BETA3(0,1)     /   -5.6265120D0/
      DATA GSS3(1)        /   14.7942080D0/
      DATA ACORE3(1)      /    3.3563860D0/
      DATA AGAUS3(1,1)    /    1.1287500D0/
      DATA AGAUS3(2,1)    /   -1.0603290D0/
      DATA BGAUS3(1,1)    /    5.0962820D0/
      DATA BGAUS3(2,1)    /    6.0037880D0/
      DATA CGAUS3(1,1)    /    1.5374650D0/
      DATA CGAUS3(2,1)    /    1.5701890D0/
      DATA EEATM3(1)      /  -13.0733210D0/
      DATA HFATM3(1)      /       52.102D0/
C
C     DATA FOR BERYLLIUM:
C
      DATA UCORE3(0,4)    /  -17.2647520D0/
      DATA UCORE3(1,4)    /  -11.3042430D0/
      DATA EXPNT3(0,4)    /    0.8774390D0/
      DATA EXPNT3(1,4)    /    1.5087550D0/
      DATA AL3(0,4)       /    0.3312330D0/
      DATA AL3(1,4)       /    0.2908996D0/
      DATA AL3(2,4)       /    0.3530008D0/
      DATA DL3(1,4)       /    1.0090531D0/
      DATA DL3(2,4)       /    0.8117586D0/
      DATA BETA3(0,4)     /   -3.9620530D0/
      DATA BETA3(1,4)     /   -2.7806840D0/
      DATA GSS3(4)        /    9.0128510D0/
      DATA GPP3(4)        /    6.0571820D0/
      DATA GSP3(4)        /    6.5761990D0/
      DATA GP23(4)        /    9.0052190D0/
      DATA HSP3(4)        /    0.5446790D0/
      DATA ACORE3(4)      /    1.5935360D0/
      DATA AGAUS3(1,4)    /    1.6315720D0/
      DATA AGAUS3(2,4)    /   -2.1109590D0/
      DATA BGAUS3(1,4)    /    2.6729620D0/
      DATA BGAUS3(2,4)    /    1.9685940D0/
      DATA CGAUS3(1,4)    /    1.7916860D0/
      DATA CGAUS3(2,4)    /    1.7558710D0/
      DATA EEATM3(4)      /  -25.5166530D0/
      DATA HFATM3(4)      /       76.960D0/
C
C     DATA FOR CARBON:
C
      DATA UCORE3(0,6)    /  -47.2703200D0/
      DATA UCORE3(1,6)    /  -36.2669180D0/
      DATA EXPNT3(0,6)    /    1.5650850D0/
      DATA EXPNT3(1,6)    /    1.8423450D0/
      DATA AL3(0,6)       /    0.4116394D0/
      DATA AL3(1,6)       /    0.5885862D0/
      DATA AL3(2,6)       /    0.7647667D0/
      DATA DL3(1,6)       /    0.8332396D0/
      DATA DL3(2,6)       /    0.6647750D0/
      DATA BETA3(0,6)     /  -11.9100150D0/
      DATA BETA3(1,6)     /   -9.8027550D0/
      DATA GSS3(6)        /   11.2007080D0/
      DATA GPP3(6)        /   10.7962920D0/
      DATA GSP3(6)        /   10.2650270D0/
      DATA GP23(6)        /    9.0425660D0/
      DATA HSP3(6)        /    2.2909800D0/
      DATA ACORE3(6)      /    2.7078070D0/
      DATA AGAUS3(1,6)    /    0.0501070D0/
      DATA AGAUS3(2,6)    /    0.0507330D0/
      DATA BGAUS3(1,6)    /    6.0031650D0/
      DATA BGAUS3(2,6)    /    6.0029790D0/
      DATA CGAUS3(1,6)    /    1.6422140D0/
      DATA CGAUS3(2,6)    /    0.8924880D0/
      DATA EEATM3(6)      / -111.2299170D0/
      DATA HFATM3(6)      /      170.890D0/
C
C     DATA FOR NITROGEN:
C
      DATA UCORE3(0,7)    /  -49.3356720D0/
      DATA UCORE3(1,7)    /  -47.5097360D0/
      DATA EXPNT3(0,7)    /    2.0280940D0/
      DATA EXPNT3(1,7)    /    2.3137280D0/
      DATA AL3(0,7)       /    0.4375151D0/
      DATA AL3(1,7)       /    0.5030995D0/
      DATA AL3(2,7)       /    0.7364933D0/
      DATA DL3(1,7)       /    0.6577006D0/
      DATA DL3(2,7)       /    0.5293383D0/
      DATA BETA3(0,7)     /  -14.0625210D0/
      DATA BETA3(1,7)     /  -20.0438480D0/
      DATA GSS3(7)        /   11.9047870D0/
      DATA GPP3(7)        /   11.7546720D0/
      DATA GSP3(7)        /    7.3485650D0/
      DATA GP23(7)        /   10.8072770D0/
      DATA HSP3(7)        /    1.1367130D0/
      DATA ACORE3(7)      /    2.8305450D0/
      DATA AGAUS3(1,7)    /    1.5016740D0/
      DATA AGAUS3(2,7)    /   -1.5057720D0/
      DATA BGAUS3(1,7)    /    5.9011480D0/
      DATA BGAUS3(2,7)    /    6.0046580D0/
      DATA CGAUS3(1,7)    /    1.7107400D0/
      DATA CGAUS3(2,7)    /    1.7161490D0/
      DATA EEATM3(7)      / -157.6137755D0/
      DATA HFATM3(7)      /      113.000D0/
C
C     DATA FOR OXYGEN:
C
      DATA UCORE3(0,8)    /  -86.9930020D0/
      DATA UCORE3(1,8)    /  -71.8795800D0/
      DATA EXPNT3(0,8)    /    3.7965440D0/
      DATA EXPNT3(1,8)    /    2.3894020D0/
      DATA AL3(0,8)       /    0.5790430D0/
      DATA AL3(1,8)       /    0.5299630D0/
      DATA AL3(2,8)       /    0.8179630D0/
      DATA DL3(1,8)       /    0.4086173D0/
      DATA DL3(2,8)       /    0.5125738D0/
      DATA BETA3(0,8)     /  -45.2026510D0/
      DATA BETA3(1,8)     /  -24.7525150D0/
      DATA GSS3(8)        /   15.7557600D0/
      DATA GPP3(8)        /   13.6540160D0/
      DATA GSP3(8)        /   10.6211600D0/
      DATA GP23(8)        /   12.4060950D0/
      DATA HSP3(8)        /    0.5938830D0/
      DATA ACORE3(8)      /    3.2171020D0/
      DATA AGAUS3(1,8)    /   -1.1311280D0/
      DATA AGAUS3(2,8)    /    1.1378910D0/
      DATA BGAUS3(1,8)    /    6.0024770D0/
      DATA BGAUS3(2,8)    /    5.9505120D0/
      DATA CGAUS3(1,8)    /    1.6073110D0/
      DATA CGAUS3(2,8)    /    1.5983950D0/
      DATA EEATM3(8)      / -289.3422065D0/
      DATA HFATM3(8)      /       59.559D0/
C
C     DATA FOR FLUORINE:
C
      DATA UCORE3(0,9)    / -110.4353030D0/
      DATA UCORE3(1,9)    / -105.6850470D0/
      DATA EXPNT3(0,9)    /    4.7085550D0/
      DATA EXPNT3(1,9)    /    2.4911780D0/
      DATA AL3(0,9)       /    0.3857650D0/
      DATA AL3(1,9)       /    0.6768503D0/
      DATA AL3(2,9)       /    0.6120047D0/
      DATA DL3(1,9)       /    0.3125302D0/
      DATA DL3(2,9)       /    0.4916328D0/
      DATA BETA3(0,9)     /  -48.4059390D0/
      DATA BETA3(1,9)     /  -27.7446600D0/
      DATA GSS3(9)        /   10.4966670D0/
      DATA GPP3(9)        /   14.8172560D0/
      DATA GSP3(9)        /   16.0736890D0/
      DATA GP23(9)        /   14.4183930D0/
      DATA HSP3(9)        /    0.7277630D0/
      DATA ACORE3(9)      /    3.3589210D0/
      DATA AGAUS3(1,9)    /   -0.0121660D0/
      DATA AGAUS3(2,9)    /   -0.0028520D0/
      DATA BGAUS3(1,9)    /    6.0235740D0/
      DATA BGAUS3(2,9)    /    6.0037170D0/
      DATA CGAUS3(1,9)    /    1.8568590D0/
      DATA CGAUS3(2,9)    /    2.6361580D0/
      DATA EEATM3(9)      / -437.5171690D0/
      DATA HFATM3(9)      /       18.890D0/
C
C     DATA FOR SODIUM-LIKE SPARKLE:
C
! 2007-06-08: change sparkle definition to Ed's Na (GM-JT)
!
!     DATA UCORE3(0,11)   /    0.0000000D0/
!     DATA UCORE3(1,11)   /    0.0000000D0/
!     DATA EXPNT3(0,11)   /    0.0000000D0/
!     DATA EXPNT3(1,11)   /    0.0000000D0/
!     DATA EXPNT3(2,11)   /    0.0000000D0/
!     DATA AL3(0,11)      /    0.5000000D0/
!     DATA AL3(1,11)      /    0.0000000D0/
!     DATA AL3(2,11)      /    0.0000000D0/
!     DATA DL3(1,11)      /    0.0000000D0/
!     DATA DL3(2,11)      /    0.0000000D0/
!     DATA BETA3(0,11)    /    0.0000000D0/
!     DATA BETA3(1,11)    /    0.0000000D0/
!     DATA GSS3(11)       /    0.0000000D0/
!     DATA GPP3(11)       /    0.0000000D0/
!     DATA GSP3(11)       /    0.0000000D0/
!     DATA GP23(11)       /    0.0000000D0/
!     DATA HSP3(11)       /    0.0000000D0/
!     DATA ACORE3(11)     /    1.6810000D0/
!     DATA AGAUS3(1,11)   /    0.0000000D0/
!     DATA AGAUS3(2,11)   /    0.0000000D0/
!     DATA BGAUS3(1,11)   /    0.0000000D0/
!     DATA BGAUS3(2,11)   /    0.0000000D0/
!     DATA CGAUS3(1,11)   /    0.0000000D0/
!     DATA CGAUS3(2,11)   /    0.0000000D0/
!     DATA EEATM3(11)     /    0.0000000D0/
!     DATA HFATM3(11)     /   25.8500000D0/

      DATA UCORE3(0,11)   /   -5.2945929D0/
      DATA UCORE3(1,11)   /   -2.4596564D0/
      DATA EXPNT3(0,11)   /    1.1375011D0/
      DATA EXPNT3(1,11)   /    1.1877433D0/
      DATA EXPNT3(2,11)   /    0.0000000D0/
      DATA BETA3(0,11)    /   -0.1510870D0/
      DATA BETA3(1,11)    /   -0.2184096D0/
      DATA GSS3(11)       /    3.9558692D0/
      DATA GPP3(11)       /    5.3363963D0/
      DATA GSP3(11)       /    7.1929109D0/
      DATA GP23(11)       /    5.0588074D0/
      DATA HSP3(11)       /    0.5687889D0/
      DATA ACORE3(11)     /    2.3677169D0/
      DATA AGAUS3(1,11)   /    0.6433655D0/
      DATA AGAUS3(2,11)   /    1.0871788D0/
      DATA BGAUS3(1,11)   /    1.5465054D0/
      DATA BGAUS3(2,11)   /    1.4529000D0/
      DATA CGAUS3(1,11)   /    0.9976699D0/
      DATA CGAUS3(2,11)   /    1.4506099D0/
      DATA DL3(1,11)      /    1.7352377D0/
      DATA DL3(2,11)      /    1.4088230D0/
      DATA AL3(0,11)      /    0.1453829D0/
      DATA AL3(1,11)      /    0.2131879D0/
      DATA AL3(2,11)      /    0.2575829D0/
      DATA EEATM3(11)     /   -5.2945929D0/
      DATA HFATM3(11)     /   25.8500000D0/

C
C     DATA FOR MAGNESIUM:
C
      DATA UCORE3(0,12)   /  -14.6236880D0/
      DATA UCORE3(1,12)   /  -14.1734600D0/
      DATA EXPNT3(0,12)   /    0.6985520D0/
      DATA EXPNT3(1,12)   /    1.4834530D0/
      DATA AL3(0,12)      /    0.2460235D0/
      DATA AL3(1,12)      /    0.2695751D0/
      DATA AL3(2,12)      /    0.2767522D0/
      DATA DL3(1,12)      /    1.1403950D0/
      DATA DL3(2,12)      /    1.1279899D0/
      DATA BETA3(0,12)    /   -2.0716910D0/
      DATA BETA3(1,12)    /   -0.5695810D0/
      DATA GSS3(12)       /    6.6943000D0/
      DATA GPP3(12)       /    6.9104460D0/
      DATA GSP3(12)       /    6.7939950D0/
      DATA GP23(12)       /    7.0908230D0/
      DATA HSP3(12)       /    0.5433000D0/
      DATA ACORE3(12)     /    1.3291470D0/
      DATA AGAUS3(1,12)   /    2.1170500D0/
      DATA AGAUS3(2,12)   /   -2.5477670D0/
      DATA BGAUS3(1,12)   /    6.0094770D0/
      DATA BGAUS3(2,12)   /    4.3953700D0/
      DATA CGAUS3(1,12)   /    2.0844060D0/
      DATA CGAUS3(2,12)   /    2.0636740D0/
      DATA EEATM3(12)     /  -22.5530760D0/
      DATA HFATM3(12)     /       35.000D0/
C
C     DATA FOR ALUMINUM:
C
      DATA UCORE3(0,13)   /  -24.8454040D0/
      DATA UCORE3(1,13)   /  -22.2641590D0/
      DATA EXPNT3(0,13)   /    1.7028880D0/
      DATA EXPNT3(1,13)   /    1.0736290D0/
      DATA EXPNT3(2,13)   /    1.0000000D0/
      DATA AL3(0,13)      /    0.2123020D0/
      DATA AL3(1,13)      /    0.6418584D0/
      DATA AL3(2,13)      /    0.2262838D0/
      DATA DL3(1,13)      /    1.2102799D0/
      DATA DL3(2,13)      /    1.5585645D0/
      DATA BETA3(0,13)    /   -0.5943010D0/
      DATA BETA3(1,13)    /   -0.9565500D0/
      DATA GSS3(13)       /    5.7767370D0/
      DATA GPP3(13)       /    6.3477900D0/
      DATA GSP3(13)       /   11.6598560D0/
      DATA GP23(13)       /    6.1210770D0/
      DATA HSP3(13)       /    4.0062450D0/
      DATA ACORE3(13)     /    1.5217030D0/
      DATA AGAUS3(1,13)   /   -0.4730900D0/
      DATA AGAUS3(2,13)   /   -0.1540510D0/
      DATA BGAUS3(1,13)   /    1.9158250D0/
      DATA BGAUS3(2,13)   /    6.0050860D0/
      DATA CGAUS3(1,13)   /    1.4517280D0/
      DATA CGAUS3(2,13)   /    2.5199970D0/
      DATA EEATM3(13)     /  -46.8647630D0/
      DATA HFATM3(13)     /       79.490D0/
C
C     DATA FOR SILICON:
C
      DATA UCORE3(0,14)   /  -26.7634830D0/
      DATA UCORE3(1,14)   /  -22.8136350D0/
      DATA EXPNT3(0,14)   /    1.6350750D0/
      DATA EXPNT3(1,14)   /    1.3130880D0/
      DATA EXPNT3(2,14)   /    1.0000000D0/
      DATA AL3(0,14)      /    0.1854905D0/
      DATA AL3(1,14)      /    0.3060715D0/
      DATA AL3(2,14)      /    0.4877432D0/
      DATA DL3(1,14)      /    1.3144550D0/
      DATA DL3(2,14)      /    1.2743396D0/
      DATA BETA3(0,14)    /   -2.8621450D0/
      DATA BETA3(1,14)    /   -3.9331480D0/
      DATA GSS3(14)       /    5.0471960D0/
      DATA GPP3(14)       /    6.7593670D0/
      DATA GSP3(14)       /    5.9490570D0/
      DATA GP23(14)       /    5.1612970D0/
      DATA HSP3(14)       /    0.9198320D0/
      DATA ACORE3(14)     /    2.1358090D0/
      DATA AGAUS3(1,14)   /   -0.3906000D0/
      DATA AGAUS3(2,14)   /    0.0572590D0/
      DATA BGAUS3(1,14)   /    6.0000540D0/
      DATA BGAUS3(2,14)   /    6.0071830D0/
      DATA CGAUS3(1,14)   /    0.6322620D0/
      DATA CGAUS3(2,14)   /    2.0199870D0/
      DATA EEATM3(14)     /  -67.7882140D0/
      DATA HFATM3(14)     /      108.390D0/
C
C     DATA FOR PHOSPHORUS:
C
      DATA UCORE3(0,15)   /  -40.4130960D0/
      DATA UCORE3(1,15)   /  -29.5930520D0/
      DATA EXPNT3(0,15)   /    2.0175630D0/
      DATA EXPNT3(1,15)   /    1.5047320D0/
      DATA EXPNT3(2,15)   /    1.0000000D0/
      DATA AL3(0,15)      /    0.2867187D0/
      DATA AL3(1,15)      /    0.4309446D0/
      DATA AL3(2,15)      /    0.3732517D0/
      DATA DL3(1,15)      /    1.0644947D0/
      DATA DL3(2,15)      /    1.1120386D0/
      DATA BETA3(0,15)    /  -12.6158790D0/
      DATA BETA3(1,15)    /   -4.1600400D0/
      DATA GSS3(15)       /    7.8016150D0/
      DATA GPP3(15)       /    6.6184780D0/
      DATA GSP3(15)       /    5.1869490D0/
      DATA GP23(15)       /    6.0620020D0/
      DATA HSP3(15)       /    1.5428090D0/
      DATA ACORE3(15)     /    1.9405340D0/
      DATA AGAUS3(1,15)   /   -0.6114210D0/
      DATA AGAUS3(2,15)   /   -0.0939350D0/
      DATA BGAUS3(1,15)   /    1.9972720D0/
      DATA BGAUS3(2,15)   /    1.9983600D0/
      DATA CGAUS3(1,15)   /    0.7946240D0/
      DATA CGAUS3(2,15)   /    1.9106770D0/
      DATA EEATM3(15)     / -117.9591740D0/
      DATA HFATM3(15)     /       75.570D0/
C
C     DATA FOR SULFUR:
C
      DATA UCORE3(0,16)   /  -49.8953710D0/
      DATA UCORE3(1,16)   /  -44.3925830D0/
      DATA EXPNT3(0,16)   /    1.8911850D0/
      DATA EXPNT3(1,16)   /    1.6589720D0/
      DATA EXPNT3(2,16)   /    1.0000000D0/
      DATA AL3(0,16)      /    0.3294622D0/
      DATA AL3(1,16)      /    0.6679118D0/
      DATA AL3(2,16)      /    0.6137472D0/
      DATA DL3(1,16)      /    1.1214313D0/
      DATA DL3(2,16)      /    1.0086488D0/
      DATA BETA3(0,16)    /   -8.8274650D0/
      DATA BETA3(1,16)    /   -8.0914150D0/
      DATA GSS3(16)       /    8.9646670D0/
      DATA GPP3(16)       /    9.9681640D0/
      DATA GSP3(16)       /    6.7859360D0/
      DATA GP23(16)       /    7.9702470D0/
      DATA HSP3(16)       /    4.0418360D0/
      DATA ACORE3(16)     /    2.2697060D0/
      DATA AGAUS3(1,16)   /   -0.3991910D0/
      DATA AGAUS3(2,16)   /   -0.0548990D0/
      DATA BGAUS3(1,16)   /    6.0006690D0/
      DATA BGAUS3(2,16)   /    6.0018450D0/
      DATA CGAUS3(1,16)   /    0.9621230D0/
      DATA CGAUS3(2,16)   /    1.5799440D0/
      DATA EEATM3(16)     / -183.4537395D0/
      DATA HFATM3(16)     /       66.400D0/
C
C     DATA FOR CHLORINE:
C
      DATA UCORE3(0,17)   / -100.6267470D0/
      DATA UCORE3(1,17)   /  -53.6143960D0/
      DATA EXPNT3(0,17)   /    2.2462100D0/
      DATA EXPNT3(1,17)   /    2.1510100D0/
      DATA EXPNT3(2,17)   /    1.0000000D0/
      DATA AL3(0,17)      /    0.5885190D0/
      DATA AL3(1,17)      /    0.6814522D0/
      DATA AL3(2,17)      /    0.3643694D0/
      DATA DL3(1,17)      /    0.9175856D0/
      DATA DL3(2,17)      /    0.7779230D0/
      DATA BETA3(0,17)    /  -27.5285600D0/
      DATA BETA3(1,17)    /  -11.5939220D0/
      DATA GSS3(17)       /   16.0136010D0/
      DATA GPP3(17)       /    7.5222150D0/
      DATA GSP3(17)       /    8.0481150D0/
      DATA GP23(17)       /    7.5041540D0/
      DATA HSP3(17)       /    3.4811530D0/
      DATA ACORE3(17)     /    2.5172960D0/
      DATA AGAUS3(1,17)   /   -0.1715910D0/
      DATA AGAUS3(2,17)   /   -0.0134580D0/
      DATA BGAUS3(1,17)   /    6.0008020D0/
      DATA BGAUS3(2,17)   /    1.9666180D0/
      DATA CGAUS3(1,17)   /    1.0875020D0/
      DATA CGAUS3(2,17)   /    2.2928910D0/
      DATA EEATM3(17)     / -315.1949480D0/
      DATA HFATM3(17)     /       28.990D0/
C
C     DATA FOR ZINC:
C
      DATA UCORE3(0,30)   /  -18.5321980D0/
      DATA UCORE3(1,30)   /  -11.0474090D0/
      DATA EXPNT3(0,30)   /    1.8199890D0/
      DATA EXPNT3(1,30)   /    1.5069220D0/
      DATA EXPNT3(2,30)   /    1.0000000D0/
      DATA AL3(0,30)      /    0.3556485D0/
      DATA AL3(1,30)      /    0.2375689D0/
      DATA AL3(2,30)      /    0.2661069D0/
      DATA DL3(1,30)      /    1.5005758D0/
      DATA DL3(2,30)      /    1.4077174D0/
      DATA BETA3(0,30)    /   -0.7155780D0/
      DATA BETA3(1,30)    /   -6.3518640D0/
      DATA GSS3(30)       /    9.6771960D0/
      DATA GPP3(30)       /    4.9801740D0/
      DATA GSP3(30)       /    7.7362040D0/
      DATA GP23(30)       /    4.6696560D0/
      DATA HSP3(30)       /    0.6004130D0/
      DATA ACORE3(30)     /    1.3501260D0/
      DATA AGAUS3(1,30)   /   -0.1112340D0/
      DATA AGAUS3(2,30)   /   -0.1323700D0/
      DATA BGAUS3(1,30)   /    6.0014780D0/
      DATA BGAUS3(2,30)   /    1.9958390D0/
      DATA CGAUS3(1,30)   /    1.5160320D0/
      DATA CGAUS3(2,30)   /    2.5196420D0/
      DATA EEATM3(30)     /  -27.3872000D0/
      DATA HFATM3(30)     /       31.170D0/
C
C     DATA FOR GALLIUM:
C
      DATA UCORE3(0,31)   /  -29.8555930D0/
      DATA UCORE3(1,31)   /  -21.8753710D0/
      DATA EXPNT3(0,31)   /    1.8470400D0/
      DATA EXPNT3(1,31)   /    0.8394110D0/
      DATA AL3(0,31)      /    0.3108620D0/
      DATA AL3(1,31)      /    0.5129360D0/
      DATA AL3(2,31)      /    0.1546208D0/
      DATA DL3(1,31)      /    0.9776692D0/
      DATA DL3(2,31)      /    2.5271534D0/
      DATA BETA3(0,31)    /   -4.9456180D0/
      DATA BETA3(1,31)    /   -0.4070530D0/
      DATA GSS3(31)       /    8.4585540D0/
      DATA GPP3(31)       /    5.0868550D0/
      DATA GSP3(31)       /    8.9256190D0/
      DATA GP23(31)       /    4.9830450D0/
      DATA HSP3(31)       /    2.0512600D0/
      DATA ACORE3(31)     /    1.6051150D0/
      DATA AGAUS3(1,31)   /   -0.5601790D0/
      DATA AGAUS3(2,31)   /   -0.2727310D0/
      DATA BGAUS3(1,31)   /    5.6232730D0/
      DATA BGAUS3(2,31)   /    1.9918430D0/
      DATA CGAUS3(1,31)   /    1.5317800D0/
      DATA CGAUS3(2,31)   /    2.1838640D0/
      DATA EEATM3(31)     /  -57.3280250D0/
      DATA HFATM3(31)     /       65.400D0/
C
C     DATA FOR GERMANIUM:
C
      DATA UCORE3(0,32)   /  -35.4671955D0/
      DATA UCORE3(1,32)   /  -31.5863583D0/
      DATA EXPNT3(0,32)   /    2.2373526D0/
      DATA EXPNT3(1,32)   /    1.5924319D0/
      DATA AL3(0,32)      /    0.1976098D0/
      DATA AL3(1,32)      /    0.3798182D0/
      DATA AL3(2,32)      /    0.3620669D0/
      DATA DL3(1,32)      /    1.1920304D0/
      DATA DL3(2,32)      /    1.3321263D0/
      DATA BETA3(0,32)    /   -5.3250024D0/
      DATA BETA3(1,32)    /   -2.2501567D0/
      DATA GSS3(32)       /    5.3769635D0/
      DATA GPP3(32)       /    7.6718647D0/
      DATA GSP3(32)       /   10.2095293D0/
      DATA GP23(32)       /    6.9242663D0/
      DATA HSP3(32)       /    1.3370204D0/
      DATA ACORE3(32)     /    1.9723370D0/
      DATA AGAUS3(1,32)   /    0.9631726D0/
      DATA AGAUS3(2,32)   /   -0.9593891D0/
      DATA BGAUS3(1,32)   /    6.0120134D0/
      DATA BGAUS3(2,32)   /    5.7491802D0/
      DATA CGAUS3(1,32)   /    2.1633655D0/
      DATA CGAUS3(2,32)   /    2.1693724D0/
      DATA EEATM3(32)     /  -84.0156006D0/
      DATA HFATM3(32)     /       89.500D0/
C
C     DATA FOR ARSENIC:
C
      DATA UCORE3(0,33)   /  -38.5074240D0/
      DATA UCORE3(1,33)   /  -35.1524150D0/
      DATA EXPNT3(0,33)   /    2.6361770D0/
      DATA EXPNT3(1,33)   /    1.7038890D0/
      DATA AL3(0,33)      /    0.3230063D0/
      DATA AL3(1,33)      /    0.5042239D0/
      DATA AL3(2,33)      /    0.2574219D0/
      DATA DL3(1,33)      /    0.9679655D0/
      DATA DL3(2,33)      /    1.2449874D0/
      DATA BETA3(0,33)    /   -8.2321650D0/
      DATA BETA3(1,33)    /   -5.0173860D0/
      DATA GSS3(33)       /    8.7890010D0/
      DATA GPP3(33)       /    8.2872500D0/
      DATA GSP3(33)       /    5.3979830D0/
      DATA GP23(33)       /    8.2103460D0/
      DATA HSP3(33)       /    1.9510340D0/
      DATA ACORE3(33)     /    1.7944770D0/
      DATA AGAUS3(1,33)   /   -0.4600950D0/
      DATA AGAUS3(2,33)   /   -0.0889960D0/
      DATA BGAUS3(1,33)   /    1.9831150D0/
      DATA BGAUS3(2,33)   /    1.9929440D0/
      DATA CGAUS3(1,33)   /    1.0867930D0/
      DATA CGAUS3(2,33)   /    2.1400580D0/
      DATA EEATM3(33)     / -122.6326140D0/
      DATA HFATM3(33)     /       72.300D0/
C
C     DATA FOR SELENIUM:
C
      DATA UCORE3(0,34)   /  -55.3781350D0/
      DATA UCORE3(1,34)   /  -49.8230760D0/
      DATA EXPNT3(0,34)   /    2.8280510D0/
      DATA EXPNT3(1,34)   /    1.7325360D0/
      DATA AL3(0,34)      /    0.2731566D0/
      DATA AL3(1,34)      /    0.7509697D0/
      DATA AL3(2,34)      /    0.5283737D0/
      DATA DL3(1,34)      /    0.8719813D0/
      DATA DL3(2,34)      /    1.2244019D0/
      DATA BETA3(0,34)    /   -6.1578220D0/
      DATA BETA3(1,34)    /   -5.4930390D0/
      DATA GSS3(34)       /    7.4325910D0/
      DATA GPP3(34)       /    9.5683260D0/
      DATA GSP3(34)       /   10.0604610D0/
      DATA GP23(34)       /    7.7242890D0/
      DATA HSP3(34)       /    4.0165580D0/
      DATA ACORE3(34)     /    3.0439570D0/
      DATA AGAUS3(1,34)   /    0.0478730D0/
      DATA AGAUS3(2,34)   /    0.1147200D0/
      DATA BGAUS3(1,34)   /    6.0074000D0/
      DATA BGAUS3(2,34)   /    6.0086720D0/
      DATA CGAUS3(1,34)   /    2.0817170D0/
      DATA CGAUS3(2,34)   /    1.5164230D0/
      DATA EEATM3(34)     / -192.7748115D0/
      DATA HFATM3(34)     /       54.300D0/
C
C     DATA FOR BROMINE:
C
      DATA UCORE3(0,35)   / -116.6193110D0/
      DATA UCORE3(1,35)   /  -74.2271290D0/
      DATA EXPNT3(0,35)   /    5.3484570D0/
      DATA EXPNT3(1,35)   /    2.1275900D0/
      DATA EXPNT3(2,35)   /    1.0000000D0/
      DATA AL3(0,35)      /    0.5859399D0/
      DATA AL3(1,35)      /    0.6755383D0/
      DATA AL3(2,35)      /    0.3823719D0/
      DATA DL3(1,35)      /    0.2759025D0/
      DATA DL3(2,35)      /    0.9970532D0/
      DATA BETA3(0,35)    /  -31.1713420D0/
      DATA BETA3(1,35)    /   -6.8140130D0/
      DATA GSS3(35)       /   15.9434250D0/
      DATA GPP3(35)       /    8.2827630D0/
      DATA GSP3(35)       /   16.0616800D0/
      DATA GP23(35)       /    7.8168490D0/
      DATA HSP3(35)       /    0.5788690D0/
      DATA ACORE3(35)     /    2.5118420D0/
      DATA AGAUS3(1,35)   /    0.9604580D0/
      DATA AGAUS3(2,35)   /   -0.9549160D0/
      DATA BGAUS3(1,35)   /    5.9765080D0/
      DATA BGAUS3(2,35)   /    5.9447030D0/
      DATA CGAUS3(1,35)   /    2.3216540D0/
      DATA CGAUS3(2,35)   /    2.3281420D0/
      DATA EEATM3(35)     / -352.5398970D0/
      DATA HFATM3(35)     /       26.740D0/
C
C     DATA FOR CADMIUM:
C
      DATA UCORE3(0,48)   /  -15.8285840D0/
      DATA UCORE3(1,48)   /    8.7497950D0/
      DATA EXPNT3(0,48)   /    1.6793510D0/
      DATA EXPNT3(1,48)   /    2.0664120D0/
      DATA AL3(0,48)      /    0.3383668D0/
      DATA AL3(1,48)      /    0.3570290D0/
      DATA AL3(2,48)      /    0.2820582D0/
      DATA DL3(1,48)      /    1.5982681D0/
      DATA DL3(2,48)      /    1.2432402D0/
      DATA BETA3(0,48)    /   -8.5819440D0/
      DATA BETA3(1,48)    /   -0.6010340D0/
      DATA GSS3(48)       /    9.2069600D0/
      DATA GPP3(48)       /    4.9481040D0/
      DATA GSP3(48)       /    8.2315390D0/
      DATA GP23(48)       /    4.6696560D0/
      DATA HSP3(48)       /    1.6562340D0/
      DATA ACORE3(48)     /    1.5253820D0/
      DATA EEATM3(48)     /  -22.4502080D0/
      DATA HFATM3(48)     /       26.720D0/
C
C     DATA FOR INDIUM:
C
      DATA UCORE3(0,49)   /  -26.1762050D0/
      DATA UCORE3(1,49)   /  -20.0058220D0/
      DATA EXPNT3(0,49)   /    2.0161160D0/
      DATA EXPNT3(1,49)   /    1.4453500D0/
      DATA AL3(0,49)      /    0.2409004D0/
      DATA AL3(1,49)      /    0.4532655D0/
      DATA AL3(2,49)      /    0.3689812D0/
      DATA DL3(1,49)      /    1.5766241D0/
      DATA DL3(2,49)      /    1.7774563D0/
      DATA BETA3(0,49)    /   -2.9933190D0/
      DATA BETA3(1,49)    /   -1.8289080D0/
      DATA GSS3(49)       /    6.5549000D0/
      DATA GPP3(49)       /    6.2992690D0/
      DATA GSP3(49)       /    8.2298730D0/
      DATA GP23(49)       /    4.9842110D0/
      DATA HSP3(49)       /    2.6314610D0/
      DATA ACORE3(49)     /    1.4183850D0/
      DATA AGAUS3(1,49)   /   -0.3431380D0/
      DATA AGAUS3(2,49)   /   -0.1095320D0/
      DATA BGAUS3(1,49)   /    1.9940340D0/
      DATA BGAUS3(2,49)   /    5.6832170D0/
      DATA CGAUS3(1,49)   /    1.6255160D0/
      DATA CGAUS3(2,49)   /    2.8670090D0/
      DATA EEATM3(49)     /  -51.9750470D0/
      DATA HFATM3(49)     /       58.000D0/
C
C     DATA FOR TIN:
C
      DATA UCORE3(0,50)   /  -34.5501920D0/
      DATA UCORE3(1,50)   /  -25.8944190D0/
      DATA EXPNT3(0,50)   /    2.3733280D0/
      DATA EXPNT3(1,50)   /    1.6382330D0/
      DATA AL3(0,50)      /    0.3744959D0/
      DATA AL3(1,50)      /    0.3218163D0/
      DATA AL3(2,50)      /    0.2832529D0/
      DATA DL3(1,50)      /    1.3120038D0/
      DATA DL3(2,50)      /    1.5681814D0/
      DATA BETA3(0,50)    /   -2.7858020D0/
      DATA BETA3(1,50)    /   -2.0059990D0/
      DATA GSS3(50)       /   10.1900330D0/
      DATA GPP3(50)       /    5.6738100D0/
      DATA GSP3(50)       /    7.2353270D0/
      DATA GP23(50)       /    5.1822140D0/
      DATA HSP3(50)       /    1.0331570D0/
      DATA ACORE3(50)     /    1.6996500D0/
      DATA AGAUS3(1,50)   /   -0.1503530D0/
      DATA AGAUS3(2,50)   /   -0.0444170D0/
      DATA BGAUS3(1,50)   /    6.0056940D0/
      DATA BGAUS3(2,50)   /    2.2573810D0/
      DATA CGAUS3(1,50)   /    1.7046420D0/
      DATA CGAUS3(2,50)   /    2.4698690D0/
      DATA EEATM3(50)     /  -78.8877790D0/
      DATA HFATM3(50)     /       72.200D0/
C
C     DATA FOR ANTIMONY:
C
      DATA UCORE3(0,51)   /  -56.4321960D0/
      DATA UCORE3(1,51)   /  -29.4349540D0/
      DATA EXPNT3(0,51)   /    2.3430390D0/
      DATA EXPNT3(1,51)   /    1.8999920D0/
      DATA AL3(0,51)      /    0.3395177D0/
      DATA AL3(1,51)      /    0.4589010D0/
      DATA AL3(2,51)      /    0.2423472D0/
      DATA DL3(1,51)      /    1.4091903D0/
      DATA DL3(2,51)      /    1.3521354D0/
      DATA BETA3(0,51)    /  -14.7942170D0/
      DATA BETA3(1,51)    /   -2.8179480D0/
      DATA GSS3(51)       /    9.2382770D0/
      DATA GPP3(51)       /    6.3500000D0/
      DATA GSP3(51)       /    5.2776800D0/
      DATA GP23(51)       /    6.2500000D0/
      DATA HSP3(51)       /    2.4244640D0/
      DATA ACORE3(51)     /    2.0343010D0/
      DATA AGAUS3(1,51)   /    3.0020280D0/
      DATA AGAUS3(2,51)   /   -0.0188920D0/
      DATA BGAUS3(1,51)   /    6.0053420D0/
      DATA BGAUS3(2,51)   /    6.0114780D0/
      DATA CGAUS3(1,51)   /    0.8530600D0/
      DATA CGAUS3(2,51)   /    2.7933110D0/
      DATA EEATM3(51)     / -148.9382890D0/
      DATA HFATM3(51)     /       63.200D0/
C
C     DATA FOR TELLURIUM:
C
      DATA UCORE3(0,52)   /  -44.9380360D0/
      DATA UCORE3(1,52)   /  -46.3140990D0/
      DATA EXPNT3(0,52)   /    4.1654920D0/
      DATA EXPNT3(1,52)   /    1.6475550D0/
      DATA AL3(0,52)      /    0.3768862D0/
      DATA AL3(1,52)      /    1.1960743D0/
      DATA AL3(2,52)      /    0.2184786D0/
      DATA DL3(1,52)      /    0.3484177D0/
      DATA DL3(2,52)      /    1.5593085D0/
      DATA BETA3(0,52)    /   -2.6651460D0/
      DATA BETA3(1,52)    /   -3.8954300D0/
      DATA GSS3(52)       /   10.2550730D0/
      DATA GPP3(52)       /    7.7775920D0/
      DATA GSP3(52)       /    8.1691450D0/
      DATA GP23(52)       /    7.7551210D0/
      DATA HSP3(52)       /    3.7724620D0/
      DATA ACORE3(52)     /    2.4850190D0/
      DATA AGAUS3(1,52)   /    0.0333910D0/
      DATA AGAUS3(2,52)   /   -1.9218670D0/
      DATA BGAUS3(1,52)   /    5.9563790D0/
      DATA BGAUS3(2,52)   /    4.9732190D0/
      DATA CGAUS3(1,52)   /    2.2775750D0/
      DATA CGAUS3(2,52)   /    0.5242430D0/
      DATA EEATM3(52)     / -168.0945925D0/
      DATA HFATM3(52)     /       47.000D0/
C
C     DATA FOR IODINE:
C
      DATA UCORE3(0,53)   /  -96.4540370D0/
      DATA UCORE3(1,53)   /  -61.0915820D0/
      DATA EXPNT3(0,53)   /    7.0010130D0/
      DATA EXPNT3(1,53)   /    2.4543540D0/
      DATA EXPNT3(2,53)   /    1.0000000D0/
      DATA AL3(0,53)      /    0.5009902D0/
      DATA AL3(1,53)      /    1.6699104D0/
      DATA AL3(2,53)      /    0.5153082D0/
      DATA DL3(1,53)      /    0.1581469D0/
      DATA DL3(2,53)      /    1.0467302D0/
      DATA BETA3(0,53)    /  -14.4942340D0/
      DATA BETA3(1,53)    /   -5.8947030D0/
      DATA GSS3(53)       /   13.6319430D0/
      DATA GPP3(53)       /    7.2883300D0/
      DATA GSP3(53)       /   14.9904060D0/
      DATA GP23(53)       /    5.9664070D0/
      DATA HSP3(53)       /    2.6300350D0/
      DATA ACORE3(53)     /    1.9901850D0/
      DATA AGAUS3(1,53)   /   -0.1314810D0/
      DATA AGAUS3(2,53)   /   -0.0368970D0/
      DATA BGAUS3(1,53)   /    5.2064170D0/
      DATA BGAUS3(2,53)   /    6.0101170D0/
      DATA CGAUS3(1,53)   /    1.7488240D0/
      DATA CGAUS3(2,53)   /    2.7103730D0/
      DATA EEATM3(53)     / -288.3160860D0/
      DATA HFATM3(53)     /       25.517D0/
C
C     DATA FOR MERCURY:
C
      DATA UCORE3(0,80)   /  -17.7622290D0/
      DATA UCORE3(1,80)   /  -18.3307510D0/
      DATA EXPNT3(0,80)   /    1.4768850D0/
      DATA EXPNT3(1,80)   /    2.4799510D0/
      DATA AL3(0,80)      /    0.2434664D0/
      DATA AL3(1,80)      /    0.4515472D0/
      DATA AL3(2,80)      /    0.2618394D0/
      DATA DL3(1,80)      /    1.2317811D0/
      DATA DL3(2,80)      /    1.2164033D0/
      DATA BETA3(0,80)    /   -3.1013650D0/
      DATA BETA3(1,80)    /   -3.4640310D0/
      DATA GSS3(80)       /    6.6247200D0/
      DATA GPP3(80)       /   14.7092830D0/
      DATA GSP3(80)       /   10.6392970D0/
      DATA GP23(80)       /   16.0007400D0/
      DATA HSP3(80)       /    2.0363110D0/
      DATA ACORE3(80)     /    1.5293770D0/
      DATA AGAUS3(1,80)   /    1.0827200D0/
      DATA AGAUS3(2,80)   /   -0.0965530D0/
      DATA BGAUS3(1,80)   /    6.4965980D0/
      DATA BGAUS3(2,80)   /    3.9262810D0/
      DATA CGAUS3(1,80)   /    1.1951460D0/
      DATA CGAUS3(2,80)   /    2.6271600D0/
      DATA EEATM3(80)     /  -28.8997380D0/
      DATA HFATM3(80)     /       14.690D0/
C
C     DATA FOR THALLIUM:
C
      DATA UCORE3(0,81)   /  -30.0531700D0/
      DATA UCORE3(1,81)   /  -26.9206370D0/
      DATA EXPNT3(0,81)   /    6.8679210D0/
      DATA EXPNT3(1,81)   /    1.9694450D0/
      DATA AL3(0,81)      /    0.3844326D0/
      DATA AL3(1,81)      /    2.5741815D0/
      DATA AL3(2,81)      /    0.2213264D0/
      DATA DL3(1,81)      /    0.0781362D0/
      DATA DL3(2,81)      /    1.5317110D0/
      DATA BETA3(0,81)    /   -1.0844950D0/
      DATA BETA3(1,81)    /   -7.9467990D0/
      DATA GSS3(81)       /   10.4604120D0/
      DATA GPP3(81)       /    4.9927850D0/
      DATA GSP3(81)       /   11.2238830D0/
      DATA GP23(81)       /    8.9627270D0/
      DATA HSP3(81)       /    2.5304060D0/
      DATA ACORE3(81)     /    1.3409510D0/
      DATA AGAUS3(1,81)   /   -1.3613990D0/
      DATA AGAUS3(2,81)   /   -0.0454010D0/
      DATA BGAUS3(1,81)   /    3.5572260D0/
      DATA BGAUS3(2,81)   /    2.3069950D0/
      DATA CGAUS3(1,81)   /    1.0928020D0/
      DATA CGAUS3(2,81)   /    2.9650290D0/
      DATA EEATM3(81)     /  -56.6492050D0/
      DATA HFATM3(81)     /       43.550D0/
C
C     DATA FOR LEAD:
C
      DATA UCORE3(0,82)   /  -30.3227560D0/
      DATA UCORE3(1,82)   /  -24.4258340D0/
      DATA EXPNT3(0,82)   /    3.1412890D0/
      DATA EXPNT3(1,82)   /    1.8924180D0/
      DATA AL3(0,82)      /    0.2576991D0/
      DATA AL3(1,82)      /    0.4527678D0/
      DATA AL3(2,82)      /    0.2150175D0/
      DATA DL3(1,82)      /    0.9866290D0/
      DATA DL3(2,82)      /    1.5940562D0/
      DATA BETA3(0,82)    /   -6.1260240D0/
      DATA BETA3(1,82)    /   -1.3954300D0/
      DATA GSS3(82)       /    7.0119920D0/
      DATA GPP3(82)       /    5.1837800D0/
      DATA GSP3(82)       /    6.7937820D0/
      DATA GP23(82)       /    5.0456510D0/
      DATA HSP3(82)       /    1.5663020D0/
      DATA ACORE3(82)     /    1.6200450D0/
      DATA AGAUS3(1,82)   /   -0.1225760D0/
      DATA AGAUS3(2,82)   /   -0.0566480D0/
      DATA BGAUS3(1,82)   /    6.0030620D0/
      DATA BGAUS3(2,82)   /    4.7437050D0/
      DATA CGAUS3(1,82)   /    1.9015970D0/
      DATA CGAUS3(2,82)   /    2.8618790D0/
      DATA EEATM3(82)     /  -73.4660775D0/
      DATA HFATM3(82)     /       46.620D0/
C
C     DATA FOR BISMUTH:
C
      DATA UCORE3(0,83)   /  -33.4959380D0/
      DATA UCORE3(1,83)   /  -35.5210260D0/
      DATA EXPNT3(0,83)   /    4.9164510D0/
      DATA EXPNT3(1,83)   /    1.9349350D0/
      DATA AL3(0,83)      /    0.1833693D0/
      DATA AL3(1,83)      /    0.6776013D0/
      DATA AL3(2,83)      /    0.2586520D0/
      DATA DL3(1,83)      /    0.2798609D0/
      DATA DL3(2,83)      /    1.5590294D0/
      DATA BETA3(0,83)    /   -5.6072830D0/
      DATA BETA3(1,83)    /   -5.8001520D0/
      DATA GSS3(83)       /    4.9894800D0/
      DATA GPP3(83)       /    8.6960070D0/
      DATA GSP3(83)       /    6.1033080D0/
      DATA GP23(83)       /    8.3354470D0/
      DATA HSP3(83)       /    0.5991220D0/
      DATA ACORE3(83)     /    1.8574310D0/
      DATA AGAUS3(1,83)   /    2.5816930D0/
      DATA AGAUS3(2,83)   /    0.0603200D0/
      DATA BGAUS3(1,83)   /    5.0940220D0/
      DATA BGAUS3(2,83)   /    6.0015380D0/
      DATA CGAUS3(1,83)   /    0.4997870D0/
      DATA CGAUS3(2,83)   /    2.4279700D0/
      DATA EEATM3(83)     / -109.2774910D0/
      DATA HFATM3(83)     /       50.100D0/

