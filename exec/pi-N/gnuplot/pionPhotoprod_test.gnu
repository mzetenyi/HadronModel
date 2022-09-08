#!/usr/bin/gnuplot
set term pdf enhanced color
#set term eps enhanced color

set style data lines

load "styles.gp"

set xlabel "plab [GeV]"
set ylabel "{/Symbol s} [{/Symbol m}b]" offset 3,0

set nologscale y

set xrange [0.1:1.1]
set yrange [0:]
set key top right
#set nokey

#set style fill transparent solid 0.3 border 
set style fill transparent pattern 6 border 

#set out "piNdilep_dsig_dM_srt1.49_rs1_newsign.pdf"
set out "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_test_pi0p.pdf"
set label 1 "(a)  {/Symbol p}^0 p" at graph 0.45, graph 0.85

plot "../experimentalData/pion_photoprod_pi0_p_tot.exp" using ($1/1000.):4:5\
        title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50",\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch" using 2:3 title "Vrancx model" ls 2,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_old" using 2:3 title "old model" ls 9,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF" using 2:3 title "Vrancx no FF." ls 2 lw 1,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF_old" using 2:3 title "old no FF." ls 9 lw 1

set out "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_test_pimp.pdf"
set label 1 "(a)  {/Symbol p}^- p" at graph 0.45, graph 0.85

plot "../experimentalData/pion_photoprod_pi0_p_tot.exp" using ($1/1000.):4:5\
        title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50",\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch" using 2:4 title "Vrancx model" ls 2,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_old" using 2:4 title "old model" ls 9,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF" using 2:4 title "Vrancx no FF." ls 2 lw 1,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF_old" using 2:4 title "old no FF." ls 9 lw 1

set out "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_test_pipn.pdf"
set label 1 "(a)  {/Symbol p}^+ n" at graph 0.45, graph 0.85

plot "../experimentalData/pion_photoprod_pi0_p_tot.exp" using ($1/1000.):4:5\
        title "experiment" with errorbars  lt 1 lw 1 lc rgbcolor "gray50",\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch" using 2:5 title "Vrancx model" ls 2,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_old" using 2:5 title "old model" ls 9,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF" using 2:5 title "Vrancx no FF." ls 2 lw 1,\
     "../results/pionPhotoprod/test/pionPhotoprod_D1232sch_noFF_old" using 2:5 title "old no FF." ls 9 lw 1
