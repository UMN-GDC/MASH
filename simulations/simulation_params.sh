#!/bin/bash

array=$(seq 1 50)
for i in $array
do 
	A=$((i/10))
	B=$((i % 10))
	echo First variable is $A, Second variable isSecond variable is  $B
done
