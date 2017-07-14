#!/bin/bash

for ii in $(seq 1 1 10);
do
	for jj in $(seq 1 1 8);
	do
		for i in $(seq 1 1 10);
		do
			for j in $(seq 1 1 8);
			do
				if [ "$ii" -eq 10 ] && [ "$i" -eq 10 ]; 	then
					./desafio ./Vsoft/VsoftSamplesDatabase/1"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/1"$i"_"$j".txt ./result.txt
					#echo "./desafio ./Vsoft/VsoftSamplesDatabase/1"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/1"$i"_"$j".txt ./result.txt"
				elif [ "$ii" -eq 10 ] && [ "$i" -ne 10 ]; 	then
					./desafio ./Vsoft/VsoftSamplesDatabase/1"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/10"$i"_"$j".txt ./result.txt
					#echo "./desafio ./Vsoft/VsoftSamplesDatabase/1"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/10"$i"_"$j".txt ./result.txt"
				elif [ "$ii" -ne 10 ] && [ "$i" -eq 10 ];	then
					./desafio ./Vsoft/VsoftSamplesDatabase/10"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/1"$i"_"$j".txt ./result.txt
					#echo "./desafio ./Vsoft/VsoftSamplesDatabase/10"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/1"$i"_"$j".txt ./result.txt"
				else
					./desafio ./Vsoft/VsoftSamplesDatabase/10"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/10"$i"_"$j".txt ./result.txt
					#echo "./desafio ./Vsoft/VsoftSamplesDatabase/10"$ii"_"$jj".txt ./Vsoft/VsoftSamplesDatabase/10"$i"_"$j".txt ./result.txt"
				fi
			done
		done
	done
done
