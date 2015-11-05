#! /bin/bash

./compile.sh 

./run/make_workspace.exe --method m135 -u --nokappa --no_syst &
./run/make_workspace.exe --method m135 -u --no_syst &
./run/make_workspace.exe --method m135 -u --use_r4 --no_syst &
./run/make_workspace.exe --method m2l  -u --lumi 1.264 --no_syst &
./run/make_workspace.exe --method m2l  --lumi 3 --no_syst &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 --sig_strength 1 &
./run/make_workspace.exe --method m1bk_nodilep --lumi 3 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 --no_syst &

wait

./run/extract_yields.exe m135_c_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m135_nc_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m135_nor4_c_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m135_nor4_nc_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m135_nor4_nokappa_c_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m135_nor4_nokappa_nc_met400_mj400_nj69_sig0_lumi0p135.root
./run/extract_yields.exe m2l_nor4_c_met400_mj400_nj69_sig0_lumi1p264.root
./run/extract_yields.exe m2l_nor4_c_met400_mj400_nj69_sig0_lumi3.root
./run/extract_yields.exe m2l_nor4_nc_met400_mj400_nj69_sig0_lumi1p264.root
./run/extract_yields.exe m2l_nor4_nc_met400_mj400_nj69_sig0_lumi3.root
./run/extract_yields.exe m1bk_nc_met400_mj400_nj69_sig1_lumi3.root
./run/extract_yields.exe m1bk_c_met400_mj400_nj69_sig1_lumi3.root
./run/extract_yields.exe m1bk_nc_met400_mj400_nj69_sig0_lumi3.root
./run/extract_yields.exe m1bk_c_met400_mj400_nj69_sig0_lumi3.root

pdflatex m135_c_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m135_nc_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m135_nor4_c_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m135_nor4_nc_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m135_nor4_nokappa_c_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m135_nor4_nokappa_nc_met400_mj400_nj69_sig0_lumi0p135_bkg_table.tex
pdflatex m2l_nor4_c_met400_mj400_nj69_sig0_lumi1p264_bkg_table.tex
pdflatex m2l_nor4_c_met400_mj400_nj69_sig0_lumi3_bkg_table.tex
pdflatex m2l_nor4_nc_met400_mj400_nj69_sig0_lumi1p264_bkg_table.tex
pdflatex m2l_nor4_nc_met400_mj400_nj69_sig0_lumi3_bkg_table.tex
pdflatex m1bk_nc_met400_mj400_nj69_sig0_lumi3_bkg_table.tex
pdflatex m1bk_nc_met400_mj400_nj69_sig1_lumi3_bkg_table.tex
pdflatex m1bk_c_met400_mj400_nj69_sig1_lumi3_bkg_table.tex
