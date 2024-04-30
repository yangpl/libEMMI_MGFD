#set terminal pngcairo size 800,800
set terminal postscript eps enhanced color font 'Helvetica,10'
set size 0.6,0.8;

set output 'survey.eps'
set key autotitle columnhead

#set grid back
set title "Acquisition geometry"
set xlabel 'X (m)'
set ylabel 'Y (m)'
set xrange [-10000:10000]
set yrange [-10000:10000]

plot "receivers.txt" using 1:2 with points ps 1 pt 8 lc rgb "black" title 'Towline', \
"sources.txt" using 1:2 with points ps 1.5 pt 7 lc rgb "black" title 'Receiver' 

