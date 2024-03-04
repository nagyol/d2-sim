gnuplot -p << eof > "$@.pdf" 2> /dev/null
set terminal pdf
set encoding utf8
set samples 15000
set key off
set ylabel 'TVD'
set xlabel 'steps'
plot for [ file in "$@" ] file w lines x1 title "$@"
eof

#set yrange [0:1]

