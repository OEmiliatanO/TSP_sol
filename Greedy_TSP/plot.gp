reset
set terminal png size 800,800
set output 'result.png'
set title "TSP result"
set autoscale
set offsets 1,1,1,1
set size ratio -1

plot "output.txt" using ($2):($3) with linespoint title ' '
