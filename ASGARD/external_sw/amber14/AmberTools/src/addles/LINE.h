c
c parameters for "rline" and geti getr getc
c line       - character of the current command line
c point(i)   - pointer for the last character of an expression i
c point(100) - pointer to end of string
c nexp       - number of expressions
c jnkf       - unit of a junk file to simplify i/o
c
      integer point(100),nexp,jnkf
      character*300 line
      common /lchr/line
      common /lint/point,nexp,jnkf
