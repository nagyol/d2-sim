#!/bin/bash
OUT=combine_`date -Isecond`.pdf
gnuplot -p << eof > ${OUT}
set terminal pdf
set samples 15000
set key off
set ylabel 'TVD'
set xlabel 'steps'
#plot for [ file in "$@" ] file w title "$@"
plot for [i=1:"$@"] "./data/".i title "Run ".i with lines
eof

open ./${OUT}


#set yrange [0:1]

