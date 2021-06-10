#!/usr/bin/gnuplot
set term pdf enhanced color
#set term eps enhanced color

set style data lines

load "styles.gp"

set xlabel "m_{inv} [GeV]"
set ylabel "d{/Symbol s}/dm_{inv} [{/Symbol m}b/GeV]" offset 3,0

set label 2 "{/Symbol \326}s = 1.49 GeV" at graph 0.4, 0.9
set label 3 "{/Symbol f}_{{/Symbol r}NR} = 0{/Symbol \260}" at graph 0.4, 0.8

set logscale y

set yrange [0.001:100]
#set key bottom left
set nokey

#set style fill transparent solid 0.3 border 
set style fill transparent pattern 6 border 

#set out "piNdilep_dsig_dM_srt1.49_rs1_newsign.pdf"
set out "../results/piN_Ndilep_dsig_dM_srt1.49_90deg.pdf"

plot "../results/piN_Ndilep_dsig_dM_srt1.49_90deg" using 1:4:5 title "all" with filledcurves fs solid 0.3 ls 8,\
     "../results/piN_Ndilep_dsig_dM_srt1.49_90deg" using 1:4 notitle ls 8 lw 1, "" using 1:5 notitle ls 8 lw 1