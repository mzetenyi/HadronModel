#set term post eps enhanced color lw 2 size 3.5, 3 25
set term post eps enhanced color 16

set nologscale y
#set data style lines
set key top right
set xlabel "E_{lab}  [GeV]"
set ylabel "{/Symbol s} [{/Symbol m}b]"
# offset 3,0

set yrange [0:]
set xrange [0.1:]
set xtics 0.1, 0.2

set output "pion_photoprod_pi0_p.eps"
set label 1 "(a)  {/Symbol p}^0 p" at graph 0.6, graph 0.7

plot \
"../experimentalData/pion_photoprod_pi0_p_tot.exp" using ($1/1000.):4:5\
 title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50", \
"../results/pionPhotoprod/pionPhotoprod_all_sch" using 2:3 title "theory" with lines lt 1 lw 2,\
"../results/pionPhotoprod/pionPhotoprod_Born" using 2:3 title "Born" with lines lt 2 lw 2

set output "pion_photoprod_pi+_n.eps"
set label 1 "(b)  {/Symbol p}^+ n" at graph 0.6, graph 0.7
set nokey

plot \
"../experimentalData/pion_photoprod_pi+_n_tot.exp" using ($1/1000.):4:5\
 title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50", \
"../results/pionPhotoprod/pionPhotoprod_all_sch" using 2:5 title "theory" with lines lt 1 lw 2,\
"../results/pionPhotoprod/pionPhotoprod_Born" using 2:5 title "Born" with lines lt 2 lw 2

set output "pion_photoprod_pi-_p.eps"
set label 1 "(c)  {/Symbol p}^- p" at graph 0.6, graph 0.7

plot \
"../experimentalData/pion_photoprod_pi-_p_tot.exp" using ($1/1000.):4:5\
 title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50", \
"../results/pionPhotoprod/pionPhotoprod_all_sch" using 2:4 title "theory" with lines lt 1 lw 2,\
"../results/pionPhotoprod/pionPhotoprod_Born" using 2:4 title "Born" with lines lt 2 lw 2
