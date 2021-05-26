set terminal png size 4000,4000
set border 0
unset xtics
unset ytics
set output 'file.png'
set style line 1 \
    linecolor rgb 'blue' \
    linetype 1 linewidth 1 \
    pointtype 7 pointsize 1
plot 'salida.txt' index 0 with linespoints ls 1 notitle
