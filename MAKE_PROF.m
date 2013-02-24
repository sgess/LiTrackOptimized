function [profile, prof_ax] = MAKE_PROF(N_drive,N_witness,Sig_drive,Sig_witness,separation,nBins)

x = randn(1,N_drive);
y = randn(1,N_witness);
d = [Sig_drive*x, Sig_witness*y+separation];
[profile,prof_ax]=hist(d,nBins);
