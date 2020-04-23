#!/bin/sh

while getopts i:o: OPT
do
  case $OPT in
    "i" ) FLG_I="TRUE" ; VALUE_I="$OPTARG" ;;
    "o" ) FLG_O="TRUE" ; VALUE_O="$OPTARG" ;;
      * ) echo "Usage: $CMDNAME [-i INPUT_FILE] [-o OUTPUT_DIRECTORY]" 1>&2
          exit 1 ;;
  esac
done

echo "Input file : $VALUE_I"
echo "Out directory : $VALUE_O"

mkdir $VALUE_O
mkdir $VALUE_O/portion

line=`wc -l $VALUE_I | awk '{print $1}'`
num=$(( $(($line / 2)) -1))

for i in `seq 0 1`
do
    python ./MakeNcRNAMatrix.py $VALUE_I ${i} $VALUE_O
done
