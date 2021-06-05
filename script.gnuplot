set terminal png size 2000,2000
set border 0
unset xtics
unset ytics
set output 'file.png'
set style line 1 \
    linecolor rgb 'blue' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1
set style line 2 \
    linecolor rgb 'red' \
    linetype 1 linewidth 1 \
    pointtype 7 pointsize 2
plot 'salida.txt' index 0 with linespoints ls 1 notitle, '' index 1 with linespoints ls 2 notitle
