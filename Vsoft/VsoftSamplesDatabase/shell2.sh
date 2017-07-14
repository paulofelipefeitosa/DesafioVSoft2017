#!/bin/bash

for i in $(seq 1 1 8);
do
	for j in $(seq 1 1 8);
	do
		echo "./desafio ./Vsoft/VsoftSamplesDatabase/101_1.txt ./Vsoft/VsoftSamplesDatabase/10"$i"_"$j".txt ./out"
	done
done
