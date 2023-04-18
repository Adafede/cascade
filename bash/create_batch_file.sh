#!/usr/bin/env bash
# -*- coding: utf-8 -*-

# A simple script to create custom mzmine batch files from template
# (use of gsed because of MacOS)

FILEPATH=$1
MODE=$2
FILENAME=$(basename "$1" .mzML)
OUTPATH=$"params/params/batch_${FILENAME}_mzmine.xml"

echo $FILEPATH &&
echo $FILENAME &&
echo $OUTPATH &&
cp params/default/batch_template_mzmine_$MODE.xml $OUTPATH &&
gsed -i "s+MYFILEPATH+$FILEPATH+g" $OUTPATH &&
gsed -i "s+MYFILENAME+$FILENAME+g" $OUTPATH &&
echo "Done"