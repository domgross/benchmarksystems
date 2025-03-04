%% HVDC parameters
Rhvdc_f=0.457;
Lhvdc_f=83.9e-3;
Chvdc_f=1.76e-6;
Chvdc_dc=62.7e-6;

rkm=3e-2;
lkm=1.05e-3;
ckm=11e-9; 
gkm=6.5e-9;

%% Hvdc connecting main and island with WT
%lineparameters
Lline=310; %(in km)
Rline=rkm*Lline;
Lline=lkm*Lline;
Cline=ckm*Lline;
Gline=gkm*Lline;

%% HVDC connecting island with WT and SC
Lline1=510; %(in km)
Rline1=rkm*Lline1;
Lline1=lkm*Lline1;
Cline1=ckm*Lline1;
Gline1=gkm*Lline1;

%% HVDC connecting island with SC and machines
Lline2=1000;%1000; %(in km)
Rline2=rkm*Lline2;
Lline2=lkm*Lline2;
Cline2=ckm*Lline2;
Gline2=gkm*Lline2;