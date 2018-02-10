#!/bin/bash

start=`date +%s`
echo "start time : ${start}"

cat learn.sh | xargs -P7 -I{} -t bash -c '{}'

end=`date +%s`

S=`expr $end - $start`
m=`expr ${S} / 60`
H=`expr ${m} / 60`
M=`expr ${m} - ${H}\*60`
S=`expr ${S} % 60`
echo "${H}:${M}:${S}"
