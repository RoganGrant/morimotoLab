#!/bin/bash
#checks md5s for all fastq files in a directory

for i in *.fastq.gz
do
  compare=`cat ${i}.md5`
  compare=${compare:0:32}
  curMD5=`md5 $i`
  if echo $curMD5 | grep -q $compare
    then echo ${i}" is valid"
  else echo ${i}" is NOT valid"
  fi
done
