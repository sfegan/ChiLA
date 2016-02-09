#!/bin/bash

logsheet=$1
utc=$2

if ( test "$utc" == "" )
then
    utc=`date -u +%Y%m%d`
fi

utclong=`echo $utc | sed -e 's/\(....\)\(..\)\(..\)/\1-\2-\3/'`

dia_dir=${HOME}/Diagnostics

${HOME}/bin/lynx -width=1000 -dump -auth=veritas:*Nizamu http://veritas.sao.arizona.edu/private/elog/VERITAS-Observer/$logsheet | sed  -n '/Message ID/,$p' > $dia_dir/logsheets/${utclong}.txt
tar zcf $dia_dir/logsheets.tar.gz -C $dia_dir logsheets
