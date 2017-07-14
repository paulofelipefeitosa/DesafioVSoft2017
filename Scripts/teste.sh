#!/bin/bash

for ii in $(seq 1 1 10);
do
	for jj in $(seq 1 1 8);
	do
		for i in $(seq 1 1 10);
		do
			for j in $(seq 1 1 8);
			do
				echo "$ii"_"$jj" "$i"_"$j"
			done
		done
	done
done
