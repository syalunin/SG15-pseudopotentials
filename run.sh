#!/bin/bash
# calculates scalar-relativistic pseudopotentials

names=( C_ONCV_PBE-1.0 B_ONCV_PBE-1.0 N_ONCV_PBE-1.0 )

rm -f pseudo/*
for name in "${names[@]}"; do
   F1=elements/$name.dat
   F2=pseudo/$name.out
   F3=pseudo/$name.upf

   src/oncvpsp.x < $F1 > $F2
   
   awk '/PSP_UPF/{flag=1; next} /END_PSP/{flag=0} flag' $F2 > $F3

   echo "$F3 written"
   rm -f $F2
done