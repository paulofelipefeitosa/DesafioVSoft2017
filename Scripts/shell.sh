#!/bin/bash

mx=4076
for i in $(seq 0 1 $mx);
do
	./getsample $i <saida> saida2
	echo "./getsample $i <saida> saida2"
	python temp2.py
	read input;
done
