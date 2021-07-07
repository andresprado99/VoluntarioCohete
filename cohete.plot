set term gif animate

set key right
set xtics; set ytics;

set output 'cohete.gif'
set size 1,1
#
set yr[-3.5:1.2]
set xr[-2:7]

Time_max=100000

do for [ii=1:Time_max:100] {plot 'mets.txt' u 2:3 every ::ii::ii w p pt 7 ps 2 lc rgb "black" t '', '' u 2:3 every ::1::ii w l lt rgb 'black' t '', '' u 4:5 every ::ii::ii w p pt 7 ps 1 lc rgb "black" t '', '' u 4:5 every ::1::ii w l lt rgb 'black' t '', '' u 6:7 every ::ii::ii  w p pt 7 ps 1 lc rgb "black"  t '', '' u 6:7 every ::1::ii w l lt rgb 'black' t '', 'cohetef.txt' u 4:5 every ::ii::ii w p pt 7 ps 1 lc rgb "yellow" t 'cohete', '' u 4:5 every ::1::ii w l lt rgb 'yellow' t '', '' u 2:3 every ::ii::ii w p pt 7 ps 2 lc rgb "blue" t '', '' u 6:7 every ::ii::ii  w p pt 7 ps 1 lc rgb "gray"  t 'Luna', '' u 6:7 every ::1::ii w l lt rgb 'gray' t ''}
