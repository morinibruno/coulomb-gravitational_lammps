set size ratio 1
set lmargin screen 0.30
set rmargin screen 0.70
set tics nomirror
set tics out
set border 0
set ticslevel 0.08
set tics font "Arial,10"
set xlabel font "Arial,12"
set ylabel font "Arial,12"
set zlabel font "Arial,12"
set zlabel offset 8.0,5.0
set key font "Arial,10"
set cblabel font "Arial,12"
set xrange[-150.2:150.2]
set yrange[-150.2:150.2]
set zrange[0.0:1.5e-5]
set cbrange[0:6e-6]
set xlabel "Distância (m)" rotate parallel
set xlabel offset 0.0,0.0 
set xtics offset 0.0,-0.3
set ylabel "Distância (m)" rotate parallel
set ytics offset 0.3,0.0
set zlabel "Potencial (Kgm/s^2)"
set zlabel offset -7.0,0.0 rotate by 90
set cblabel "Potencial (Kgm/s^2)"
set cblabel offset 3.0,0.0 rotate by 270
set palette defined (0 "#000004", 1 "#1f0c48", 2 "#550f6d", 3 "#88226a", 4 "#a83655", 5 "#e35933", 6 "#f9950a", 7 "#f8c932", 8 "#fcffa4")
splot "s30.0-30.0-30.0p34191mapdata.csv" using 1:2:($3*(-1)) with pm3d at st notitle