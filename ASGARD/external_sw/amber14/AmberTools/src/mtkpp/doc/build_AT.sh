#!/bin/sh

echo " "
echo " ... BUILDING MTK++ Documentation ... "
echo " "
name=MTKpp_AT

# preclean
rm -f *.dvi *.log *.lof *.lot *.aux *.toc *.out *.bbl *.blg *.ilg *.ind *.ist *.idx *.loa *.not

echo " "
echo " ... RUNNING LATEX/BIBTEX ... "
echo " "
if which latex >/dev/null; then
  latex $name

  if which bibtex >/dev/null; then
    bibtex $name
    bibtex $name
  else
    echo bibtex does not exist
    exit 3
  fi

  latex $name
  latex $name
else
  echo latex does not exist
  exit 2
fi

echo " "
echo " ... RUNNING DVIPS ... "
echo " "
if which dvips >/dev/null; then
  dvips -o $name.ps $name.dvi
else
  echo dvips does not exist
  exit 4
fi

echo " "
echo " ... PS2PDF ... "
echo " "
if which ps2pdf >/dev/null; then
  ps2pdf $name.ps
else
  echo ps2pdf does not exist
  exit 5
fi

# postclean
rm -f *.dvi *.log *.lof *.lot *.aux *.toc *.out *.bbl *.blg *.ilg *.ind *.ist *.idx *.loa *.not *.ps
#
#
#
