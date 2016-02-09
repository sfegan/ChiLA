#!/bin/bash

laser=$1
logsheet=$2
utc=$3

if ( test "$utc" == "" )
then
    utc=`date -u +%Y%m%d`
    vbf_dir=${HOME}/Scratch/VBF/current/
else
    if ( test $utc -lt 20061224 )
    then
        vbf_dir=/net/gamma3/srv/raid1/data/d${utc}
    elif ( test $utc -lt 20070615 )
    then
        vbf_dir=/net/gamma3/srv/raid2/data/d${utc}
    elif ( test $utc -lt 20071111 )
    then
        vbf_dir=/net/gamma3/srv/raid3/data/d${utc}
    elif ( test $utc -lt 20080308 )
    then
        vbf_dir=/net/gamma3/srv/raid6/data/d${utc}
    else
        vbf_dir=/net/gamma3/srv/raid7/data/d${utc}
    fi
fi

if ( test "$4" != "" )
then
    vbf_dir=$4
fi

if ( test "$5" != "" )
then
  overwrite=$5
else
  overwrite=1
fi

utclong=`echo $utc | sed -e 's/\(....\)\(..\)\(..\)/\1-\2-\3/'`

red_dir=${HOME}/Scratch/Reduced
#dia_dir=${HOME}/Scratch/Diagnostics
dia_dir=/net/gamma3/srv/raid1/Diagnostics
code=${HOME}/ChiLA/VSScripts
script=${HOME}/Scratch/bin

mkdir $red_dir/d$utc
( cd $red_dir/d$utc; $code/analyze_directory.sh "$vbf_dir" . "$laser" "$overwrite" )
#$code/analyze_directory.sh $vbf_dir $red_dir/d$utc $laser )
#${script}/current_diagnostics_helper.sh $red_dir/d$utc $dia_dir/d$utc
${script}/current_diagnostics_helper_one.sh $red_dir/d$utc $dia_dir/d$utc $overwrite
#ssh gamma4 $code/analyze_directory.sh $vbf_dir $red_dir/d$utc $laser
#ssh gamma3 ${script}/current_diagnostics_helper.sh $red_dir/d$utc $dia_dir/d$utc

cat $dia_dir/d????????/diagnostics.dat > $dia_dir/diagnostics.dat

if ( test "$logsheet" != "" )
then
    echo $utc $laser >> $dia_dir/laser.dat
    ${HOME}/bin/lynx -width=1000 -dump -auth=veritas:*Nizamu http://veritas.sao.arizona.edu/private/elog/VERITAS-Observer/$logsheet | sed  -n '/Message ID/,$p' > $dia_dir/logsheets/${utclong}.txt
tar zcf $dia_dir/logsheets.tar.gz -C $dia_dir logsheets
fi
