reset
set terminal png size 800,800
set output 'result.png'
set title "TSP result"

plot "output.txt" using ($2):($3) with linespoint title ' '
