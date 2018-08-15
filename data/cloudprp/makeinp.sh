#!/bin/sh

rm -f doall.sh
for f in `ls -1 ./input/*.wvl | gawk '{n=split($1,a,"/");print a[n]}'`; do
  ff="./input/"$f
  prefix=`gawk 'BEGIN{split("'$f'",a,".wvl");print a[1]}'`
  echo ... preparing $prefix
  gawk -f ./makeinp.awk -vprefix=$prefix $ff
  echo sh ./tmp/$prefix.sh >> doall.sh
done

