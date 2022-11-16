      subroutine se_rdelem(symbol2, ierror)
      implicit none
#include "sebomd.dim"
#include "se_fulsys.h"

! convert AMBER atom types into atomic number for SEBOMD
!
! input: symbol2 the list of atom types (natoms length)
! output: iatnum the atomic number
!         ierror a flag (= 0 -> no error)

      character*4 symbol2(*)
      integer ierror

      character*4 element
      integer i

      ! search element for each atom type, then convert to atomic number
      do i = 1, natoms
        call se_rdelem_get_element(symbol2(i), element, ierror)
        if (ierror.ne.0) return
        call se_rdelem_get_atomic_number(element, iatnum(i), ierror)
        if (ierror.ne.0) return
      end do
      return
      end
!------------------------------------------------------------------------
      subroutine se_rdelem_get_element(symbol, element, ierror)
      implicit none

      character*4 symbol
      character*4 element
      integer ierror

      integer atom_type_size
      parameter (atom_type_size=214)

      character*4 type_conversion_table(2,atom_type_size)
      data type_conversion_table /
     . 'H   ', 'H   ', !   1
     . 'H0  ', 'H   ', !   2
     . 'H1  ', 'H   ', !   3
     . 'H2  ', 'H   ', !   4
     . 'H3  ', 'H   ', !   5
     . 'H4  ', 'H   ', !   6
     . 'H5  ', 'H   ', !   7
     . 'HA  ', 'H   ', !   8
     . 'HC  ', 'H   ', !   9
     . 'HD  ', 'H   ', !  10
     . 'HK  ', 'H   ', !  11
     . 'HN  ', 'H   ', !  12
     . 'HO  ', 'H   ', !  13
     . 'HP  ', 'H   ', !  14
     . 'HS  ', 'H   ', !  15
     . 'HT  ', 'H   ', !  16
     . 'HW  ', 'H   ', !  17
     . 'HZ  ', 'H   ', !  18
     . 'Ha  ', 'H   ', !  19
     . 'Hc  ', 'H   ', !  20
     . 'Ho  ', 'H   ', !  21
     . 'Hp  ', 'H   ', !  22
     . 'h1  ', 'H   ', !  23
     . 'h2  ', 'H   ', !  24
     . 'h3  ', 'H   ', !  25
     . 'h4  ', 'H   ', !  26
     . 'h5  ', 'H   ', !  27
     . 'hA  ', 'H   ', !  28
     . 'hB  ', 'H   ', !  29
     . 'hE  ', 'H   ', !  30
     . 'hL  ', 'H   ', !  31
     . 'hN  ', 'H   ', !  32
     . 'hO  ', 'H   ', !  33
     . 'hR  ', 'H   ', !  34
     . 'hS  ', 'H   ', !  35
     . 'hX  ', 'H   ', !  36
     . 'ha  ', 'H   ', !  37
     . 'hc  ', 'H   ', !  38
     . 'hn  ', 'H   ', !  39
     . 'ho  ', 'H   ', !  40
     . 'hp  ', 'H   ', !  41
     . 'hs  ', 'H   ', !  42
     . 'hw  ', 'H   ', !  43
     . 'hx  ', 'H   ', !  44
     . '2C  ', 'C   ', !  45
     . '3C  ', 'C   ', !  46
     . 'C   ', 'C   ', !  47
     . 'C*  ', 'C   ', !  48
     . 'C1  ', 'C   ', !  49
     . 'C2  ', 'C   ', !  50
     . 'C3  ', 'C   ', !  51
     . 'C4  ', 'C   ', !  52
     . 'C5  ', 'C   ', !  53
     . 'C6  ', 'C   ', !  54
     . 'C8  ', 'C   ', !  55
     . 'CA  ', 'C   ', !  56
     . 'CB  ', 'C   ', !  57
     . 'CC  ', 'C   ', !  58
     . 'CD  ', 'C   ', !  59
     . 'CG  ', 'C   ', !  60
     . 'CI  ', 'C   ', !  61
     . 'CJ  ', 'C   ', !  62
     . 'CK  ', 'C   ', !  63
     . 'CM  ', 'C   ', !  64
     . 'CN  ', 'C   ', !  65
     . 'CO  ', 'C   ', !  66
     . 'CP  ', 'C   ', !  67
     . 'CQ  ', 'C   ', !  68
     . 'CR  ', 'C   ', !  69
     . 'CS  ', 'C   ', !  70
     . 'CT  ', 'C   ', !  71
     . 'CV  ', 'C   ', !  72
     . 'CW  ', 'C   ', !  73
     . 'CX  ', 'C   ', !  74
     . 'CY  ', 'C   ', !  75
     . 'CZ  ', 'C   ', !  76
     . 'Cg  ', 'C   ', !  77
     . 'Cj  ', 'C   ', !  78
     . 'Ck  ', 'C   ', !  79
     . 'Cp  ', 'C   ', !  80
     . 'Cy  ', 'C   ', !  81
     . 'TG  ', 'C   ', !  82
     . 'c   ', 'C   ', !  83
     . 'c1  ', 'C   ', !  84
     . 'c2  ', 'C   ', !  85
     . 'c3  ', 'C   ', !  86
     . 'cA  ', 'C   ', !  87
     . 'cB  ', 'C   ', !  88
     . 'cC  ', 'C   ', !  89
     . 'cD  ', 'C   ', !  90
     . 'cP  ', 'C   ', !  91
     . 'cR  ', 'C   ', !  92
     . 'ca  ', 'C   ', !  93
     . 'cc  ', 'C   ', !  94
     . 'cd  ', 'C   ', !  95
     . 'ce  ', 'C   ', !  96
     . 'cf  ', 'C   ', !  97
     . 'cg  ', 'C   ', !  98
     . 'ch  ', 'C   ', !  99
     . 'cp  ', 'C   ', ! 100
     . 'cq  ', 'C   ', ! 101
     . 'cu  ', 'C   ', ! 102
     . 'cv  ', 'C   ', ! 103
     . 'cx  ', 'C   ', ! 104
     . 'cy  ', 'C   ', ! 105
     . 'cz  ', 'C   ', ! 106
     . 'N   ', 'N   ', ! 107
     . 'N*  ', 'N   ', ! 108
     . 'N2  ', 'N   ', ! 109
     . 'N3  ', 'N   ', ! 110
     . 'NA  ', 'N   ', ! 111
     . 'NB  ', 'N   ', ! 112
     . 'NC  ', 'N   ', ! 113
     . 'ND  ', 'N   ', ! 114
     . 'NL  ', 'N   ', ! 115
     . 'NO  ', 'N   ', ! 116
     . 'NP  ', 'N   ', ! 117
     . 'NT  ', 'N   ', ! 118
     . 'NX  ', 'N   ', ! 119
     . 'NY  ', 'N   ', ! 120
     . 'NZ  ', 'N   ', ! 121
     . 'Ng  ', 'N   ', ! 122
     . 'TN  ', 'N   ', ! 123
     . 'n   ', 'N   ', ! 124
     . 'n1  ', 'N   ', ! 125
     . 'n2  ', 'N   ', ! 126
     . 'n3  ', 'N   ', ! 127
     . 'n4  ', 'N   ', ! 128
     . 'nA  ', 'N   ', ! 129
     . 'na  ', 'N   ', ! 130
     . 'nb  ', 'N   ', ! 131
     . 'nc  ', 'N   ', ! 132
     . 'nd  ', 'N   ', ! 133
     . 'ne  ', 'N   ', ! 134
     . 'nf  ', 'N   ', ! 135
     . 'nh  ', 'N   ', ! 136
     . 'no  ', 'N   ', ! 137
     . 'O   ', 'O   ', ! 138
     . 'O2  ', 'O   ', ! 139
     . 'O3  ', 'O   ', ! 140
     . 'O4  ', 'O   ', ! 141
     . 'OA  ', 'O   ', ! 142
     . 'OD  ', 'O   ', ! 143
     . 'OH  ', 'O   ', ! 144
     . 'OM  ', 'O   ', ! 145
     . 'OP  ', 'O   ', ! 146
     . 'OQ  ', 'O   ', ! 147
     . 'OR  ', 'O   ', ! 148
     . 'OS  ', 'O   ', ! 149
     . 'OT  ', 'O   ', ! 150
     . 'OV  ', 'O   ', ! 151
     . 'OW  ', 'O   ', ! 152
     . 'OX  ', 'O   ', ! 153
     . 'OZ  ', 'O   ', ! 154
     . 'Oh  ', 'O   ', ! 155
     . 'Os  ', 'O   ', ! 156
     . 'Oy  ', 'O   ', ! 157
     . 'W   ', 'O   ', ! 158
     . 'o   ', 'O   ', ! 159
     . 'oC  ', 'O   ', ! 160
     . 'oH  ', 'O   ', ! 161
     . 'oO  ', 'O   ', ! 162
     . 'oP  ', 'O   ', ! 163
     . 'oR  ', 'O   ', ! 164
     . 'oS  ', 'O   ', ! 165
     . 'oT  ', 'O   ', ! 166
     . 'oh  ', 'O   ', ! 167
     . 'os  ', 'O   ', ! 168
     . 'ow  ', 'O   ', ! 169
     . 'F   ', 'F   ', ! 170
     . 'F-  ', 'F   ', ! 171
     . 'f   ', 'F   ', ! 172
     . 'P   ', 'P   ', ! 173
     . 'p2  ', 'P   ', ! 174
     . 'p3  ', 'P   ', ! 175
     . 'p4  ', 'P   ', ! 176
     . 'p5  ', 'P   ', ! 177
     . 'pA  ', 'P   ', ! 178
     . 'pb  ', 'P   ', ! 179
     . 'pc  ', 'P   ', ! 180
     . 'pd  ', 'P   ', ! 181
     . 'pe  ', 'P   ', ! 182
     . 'pf  ', 'P   ', ! 183
     . 'px  ', 'P   ', ! 184
     . 'py  ', 'P   ', ! 185
     . 'S   ', 'S   ', ! 186
     . 'S*  ', 'S   ', ! 187
     . 'S4  ', 'S   ', ! 188
     . 'SF  ', 'S   ', ! 189
     . 'SH  ', 'S   ', ! 190
     . 'SO  ', 'S   ', ! 191
     . 'SS  ', 'S   ', ! 192
     . 'Sm  ', 'S   ', ! 193
     . 's   ', 'S   ', ! 194
     . 's2  ', 'S   ', ! 195
     . 's4  ', 'S   ', ! 196
     . 's6  ', 'S   ', ! 197
     . 'sh  ', 'S   ', ! 198
     . 'ss  ', 'S   ', ! 199
     . 'sx  ', 'S   ', ! 200
     . 'sy  ', 'S   ', ! 201
     . 'IP  ', 'Na  ', ! 202
     . 'CL  ', 'Cl  ', ! 203
     . 'Cl  ', 'Cl  ', ! 204
     . 'Cl- ', 'Cl  ', ! 205
     . 'IM  ', 'Cl  ', ! 206
     . 'cl  ', 'Cl  ', ! 207
     . 'BR  ', 'Br  ', ! 208
     . 'Br  ', 'Br  ', ! 209
     . 'Br- ', 'Br  ', ! 210
     . 'br  ', 'Br  ', ! 211
     . 'I   ', 'I   ', ! 212
     . 'I-  ', 'I   ', ! 213
     . 'i   ', 'I   '  ! 214
     ./

      integer j

      do j = 1, atom_type_size
        if (symbol.eq.type_conversion_table(1,j)) then
          element = type_conversion_table(2,j)
          return
        end if
      end do
      write(6,*) 'ERROR in rdelem_get_element: symbol "',
     . symbol, '" unknown'
      ierror = 1
      return
      end
!------------------------------------------------------------------------
      subroutine se_rdelem_get_atomic_number(element, 
     &                                       atomic_number, 
     &                                       ierr)
      implicit none
      integer atomic_number
      integer ierr
      character*4 element

      integer number_of_known_element
      parameter (number_of_known_element = 11)
      character*4 element_to_atomic_number1(number_of_known_element)
      integer     element_to_atomic_number2(number_of_known_element)
      data element_to_atomic_number1/
     . 'H   ',
     . 'C   ',
     . 'N   ',
     . 'O   ',
     . 'F   ',
     . 'Na  ',
     . 'P   ',
     . 'S   ',
     . 'Cl  ',
     . 'Br  ',
     . 'I   '/
      data element_to_atomic_number2/
     . 1,
     . 6,
     . 7,
     . 8,
     . 9,
     . 11,
     . 15,
     . 16,
     . 17,
     . 35,
     . 53/
  
      integer j

      do j = 1, number_of_known_element
         if (element.eq.element_to_atomic_number1(j)) then
            atomic_number = element_to_atomic_number2(j)
            return
         end if
      end do
      write(6,*) 'ERROR in rdelem_get_atomic_number: element "',
     .           element,'" unknown'
      ierr = 1
      return
      end
!------------------------------------------------------------------------
