#!/usr/bin/gnuplot
set term pdf enhanced color
#set term eps enhanced color

set style data lines

load "styles.gp"

set xlabel "m_{inv} [GeV]"
set ylabel "d{/Symbol s}/dm_{inv} [mb/GeV]" offset 3,0

set label 2 "{/Symbol \326}s = 1.49 GeV" at graph 0.4, 0.9
set label 3 "{/Symbol f}_{{/Symbol r}NR} = 90{/Symbol \260}" at graph 0.4, 0.8

set nologscale y

set yrange [0:10]
set xrange [0.2:0.6]

#set key bottom left
set nokey

#set style fill transparent solid 0.3 border 
set style fill transparent pattern 6 border 

#set out "piNdilep_dsig_dM_srt1.49_rs1_newsign.pdf"
set out "../results/piN_Npipi_dsig_dM_srt1.49_90deg.pdf"

plot "../results/piN_Npipi_dsig_dM_srt1.49_90deg" using 1:4:5 title "all" with filledcurves fs solid 0.3 ls 8,\
     "../results/piN_Npipi_dsig_dM_srt1.49_90deg" using 1:4 notitle ls 8 lw 1, "" using 1:5 notitle ls 8 lw 1

set out "../results/piN_Npipi_dsig_dM_srt1.49_90deg_contribs.pdf"
set key bottom left

plot "../results/piN_Npipi_dsig_dM_srt1.49_Born" using 1:2 title "Born" ls 2 lw 2,\
     "../results/piN_Npipi_dsig_dM_srt1.49_N1520" using 1:2 title "N(1520)" ls 9 lw 2,\
     "../results/piN_Npipi_dsig_dM_srt1.49_N1535" using 1:2 title "N(1535)" ls 4 lw 2,\
     "../results/piN_Npipi_dsig_dM_srt1.49_N1440" using 1:2 title "N(1440)" ls 1 lw 2,\
     "../results/piN_Npipi_dsig_dM_srt1.49_all_90deg" using 1:2 title "all" ls 8 lw 2
