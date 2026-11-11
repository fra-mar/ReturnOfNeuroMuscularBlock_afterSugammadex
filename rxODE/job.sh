#!/bin/bash


rm -rf simulations
rm -rf outputs
rm -rf plots

timeStart=$(date +%Y-%m-%d\ %H:%M:%S)

echo -------------------
date
echo Run simulatorMain_251005
n=$(cat scenarios.csv | wc -l)

for i in $(seq 1 $(($n - 1)))
do
  echo ------------------------------
	echo $i
	Rscript simulatorMain_251005.R $i 300
	date
done

echo -------------------
date 
echo Run tofAnalyzer.R
Rscript tofAnalyzer.R

echo -------------------
date 
echo Run summarizer
Rscript summarizer.R

echo ------------------
date
echo Run plotter
Rscript plotter.R
mkdir plots
cp ./simulations/scenario_*/*.png ./plots
echo ------------------
echo START $timeStart
echo END 
date 
