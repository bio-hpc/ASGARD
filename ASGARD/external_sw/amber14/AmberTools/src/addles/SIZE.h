c
c define the lengths
c
c
c number of atoms, number atom types, number of residues
c
        integer maxnatom,maxntypes,maxnres
        parameter (maxnatom=50000,maxntypes=100,maxnres=100000)
c
c number of bond types, angle types and dihedral types
c
        integer maxbndt,maxangt,maxdiht
        parameter (maxbndt=1000,maxangt=1000,maxdiht=1000)
c
c number of bonds, angles and dihedrals of each type
c
        integer maxbnd,maxang,maxdih
        parameter (maxbnd=100000,maxang=800000,maxdih=15000)
c
c max natyp (is that used  by amber?), number of exclusions, hbond pair types
c
        integer maxnatyp,maxnext,maxnphb
        parameter (maxnatyp=100,maxnext=50000000,maxnphb=1000)
c
c max atoms per solvent molecule
c
        integer maxnspm
        parameter (maxnspm=200000)
c
c max number of copies you can make
c

        integer maxcopy
        parameter (maxcopy=80)
