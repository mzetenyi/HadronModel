#!/usr/bin/gnuplot
set term pdf enhanced color
#set term eps enhanced color

set style data lines

load "styles.gp"

set xlabel "sqrt(s) [GeV]"
set ylabel "{/Symbol s} [{/Symbol m}b]" offset 3,0

set nologscale y

#set yrange [0.001:100]
#set key bottom left
set nokey

#set style fill transparent solid 0.3 border 
set style fill transparent pattern 6 border 

#set out "piNdilep_dsig_dM_srt1.49_rs1_newsign.pdf"
set out "../results/pionPhotoprod.pdf"

plot "../results/pionPhotoprod_all" using 1:4:5 title "all" with filledcurves fs solid 0.3 ls 8,\
     "../results/pionPhotoprod_all" using 1:4 notitle ls 8 lw 1, "" using 1:5 notitle ls 8 lw 1