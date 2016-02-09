#!/bin/bash

red_dir=$1
dia_dir=$2

if ( test "$3" != "" )
then
  overwrite=$3
else
  overwrite=1
fi

mkdir $dia_dir
cd $dia_dir

if( test "$overwrite" != 0 )
then
  rm -f diagnostics.dat
fi

for f in $red_dir/x*_s2.h5
do
  psfile=`basename ${f/_s2.h5/_diagnostics.ps}`
  pdffile=${psfile/.ps/.pdf}
  if( test "$overwrite" != 0 -o ! -f ${pdffile} )
  then
    /veritas/matlab/bin/matlab -nodisplay -memmgr compact -nojvm -nodesktop -r "try; one_diagnostics('$f'); catch; msg=lasterror; disp(msg.message); disp(msg.identifier); end; exit;"
#    d=`basename $f`
#    ${HOME}/ChiLA/VSDataReduction/extract_diagnostics -o $d $f
#    /veritas/matlab/bin/matlab -memmgr compact -nojvm -nodesktop -r "try; one_diagnostics('$d'); catch; msg=lasterror; disp(msg.message); disp(msg.identifier); end; exit;"
#    rm -f $d
    if ( test -f $psfile )
    then
      ps2pdf $psfile
      rm -f $psfile
    fi
  fi
done
