#!/bin/bash
# by Denis Korablev for all Echidna folks

file_out=/home/production/run_validation_out.txt
file_bkp=/home/production/run_validation_out.txt.bkp
file_tmp=./file_tmp

cp $file_out $file_bkp
last=`wc -l $file_out|awk  '{ print $1 }'`
head -2 $file_out > $file_tmp
sed -n "3,"$last" p" $file_out|sort -n|uniq >> $file_tmp
cp $file_tmp $file_out
rm -rf $file_tmp

