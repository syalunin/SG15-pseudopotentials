#!/bin/bash
# calculates scalar-relativistic pseudopotentials

rm -f pseudo_all/*
for file in elements/*.dat; do
   name=$(basename "$file" .dat)
   F1=elements/$name.dat
   F2=pseudo_all/$name.out
   F3=pseudo_all/$name.upf

   src/oncvpsp.x < $F1 > $F2

   awk '/PSP_UPF/{flag=1; next} /END_PSP/{flag=0} flag' $F2 > $F3

   echo "$F3 written"
   rm -f $F2
done