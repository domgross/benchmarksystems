%a single PV cell modules
Ncell=60; %number of modules in one cell
Vpv_oc=40.79; % open circuit voltage of a sinlge unit
Ipv_sc=10.06; % short circuit current of a single unit 

Vpv_mpp=32.73; %maximal power point voltage per module
Ipv_mpp=9.32; %maximal power point current per module
Ppv_mpp=Vpv_mpp*Ipv_mpp;

%Diode and circuit parameters (single PV)
k_Boltz=1.3806e-23;
qe_charge=1.6022e-19;

Diode_ideality=0.97418;
Rpv_sh=118.4621; %Shunt resistance [Ohm]
Rpv_s=0.3749; %Series resistance [Ohm]

TK=273.15+25; %change 25 to operating temperature of the PV in celsius
G=1000; %Irradiation nominal point

I0=1.5579e-11; % current saturation
IL=10.0918; %light generated current for nominal conditions (that's ideal current source)
             %IL IL=(Isc+ki*(T-298.15))G/1000

VT=k_Boltz*TK*Diode_ideality/qe_charge*Ncell; %voltage needed for diod V-I curve

%% main land PV
Npar=1200; %number of parallel strings
Nser=90; %number of series connected modules per string

I0_total=Npar*I0;
IL_total=Npar*IL;
VT_total=VT*Nser;
Rpvs_total=Rpv_s*Nser/Npar;

Vpv_mppmax=Nser*Vpv_mpp; 
Ipv_mppmax=Npar*Ipv_mpp;

Vpv_ocmax=Nser*Vpv_oc;
Ipv_scmax=Npar*Ipv_sc;
Ppv_mppmax=Vpv_mppmax*Ipv_mppmax;

Rss=Vpv_ocmax/(IL_total-I0_total*(exp(Vpv_ocmax/VT_total) -1));
Rpvsh_total=Rss; %Vpv_ocmax/(IL_total-Ipv_scmax);%Rpv_sh*Nser/Npar; %check this one
%Rpvsh_total=Rpv_sh*Nser/Npar;
%Rpv_dc=(Vpv_mppmax/(0.05*(S_b)/Vpv_base));


%% island PV1
Npar1=5000; %number of parallel s trings
Nser1=90; %number of series connected modules per string

I0_total1=Npar1*I0;
IL_total1=Npar1*IL;
VT_total1=VT*Nser1;
Rpvs_total1=Rpv_s*Nser1/Npar1;

Vpv_mppmax1=Nser1*Vpv_mpp; 
Ipv_mppmax1=Npar1*Ipv_mpp;

Vpv_ocmax1=Nser1*Vpv_oc;
Ipv_scmax1=Npar1*Ipv_sc;
Ppv_mppmax1=Vpv_mppmax1*Ipv_mppmax1;

Rss1=Vpv_ocmax1/(IL_total1-I0_total1*(exp(Vpv_ocmax1/VT_total1) -1));
Rpvsh_total1=Rss1;  %Vpv_ocmax1/(IL_total1-Ipv_scmax1);%Rpv_sh*Nser/Npar; %check this one
%Rpvsh_total1=Rpv_sh*Nser1/Npar1;


%Rpv_dc1=(Vpv_base/(0.05*(S_b)/Vpv_base));

%% island PV2
Npar2=3000; %number of parallel s trings
Nser2=100; %number of series connected modules per string

I0_total2=Npar2*I0;
IL_total2=Npar2*IL;
VT_total2=VT*Nser2;
Rpvs_total2=Rpv_s*Nser2/Npar2;


% Parameters for PVs connected in series, and then those series in parallel
Vpv_mppmax2=Nser2*Vpv_mpp; 
Ipv_mppmax2=Npar2*Ipv_mpp;
Ppv_mppmax2=Vpv_mppmax2*Ipv_mppmax2;

Vpv_ocmax2=Nser2*Vpv_oc;
Ipv_scmax2=Npar2*Ipv_sc;

Rss2=Vpv_ocmax2/(IL_total2-I0_total2*(exp(Vpv_ocmax2/VT_total2) -1));
Rpvsh_total2=Rss2; %Vpv_ocmax2/(IL_total2-Ipv_scmax2);%Rpv_sh*Nser/Npar; %check this one
%Rpvsh_total2=Rpv_sh*Nser2/Npar2;
%Rpv_dc1=(Vpv_base/(0.05*(S_b)/Vpv_base));

%% MPPT PV3 in the area 3

Npar3=5500;
Nser3=95;

%Diode and circuit parameters for overall PV
I0_total3=Npar3*I0;
IL_total3=Npar3*IL;
VT_total3=VT*Nser3;
Rpvs_total3=Rpv_s*Nser3/Npar3;


% Parameters for PVs connected in series, and then those series in parallel
Vpv_mppmax3=Nser3*Vpv_mpp; 
Ipv_mppmax3=Npar3*Ipv_mpp;
Ppv_mppmax3=Vpv_mppmax3*Ipv_mppmax3;

Vpv_ocmax3=Nser3*Vpv_oc;
Ipv_scmax3=Npar3*Ipv_sc;

Rss3=Vpv_ocmax3/(IL_total3-I0_total3*(exp(Vpv_ocmax3/VT_total3) -1));
Rpvsh_total3=Rss3; %Vpv_ocmax2/(IL_total2-Ipv_scmax2);%Rpv_sh*Nser/Npar; %check this one
%Rpvsh_total2=Rpv_sh*Nser2/Npar2;

%% MPPT PV4 in the area 3

Npar4=5700;
Nser4=100;

%Diode and circuit parameters for overall PV
I0_total4=Npar4*I0;
IL_total4=Npar4*IL;
VT_total4=VT*Nser4;
Rpvs_total4=Rpv_s*Nser4/Npar4;


% Parameters for PVs connected in series, and then those series in parallel
Vpv_mppmax4=Nser4*Vpv_mpp; 
Ipv_mppmax4=Npar4*Ipv_mpp;
Ppv_mppmax4=Vpv_mppmax4*Ipv_mppmax4;

Vpv_ocmax4=Nser4*Vpv_oc;
Ipv_scmax4=Npar4*Ipv_sc;

Rss4=Vpv_ocmax4/(IL_total4-I0_total4*(exp(Vpv_ocmax4/VT_total4) -1));
Rpvsh_total4=Rss4; %Vpv_ocmax2/(IL_total2-Ipv_scmax2);%Rpv_sh*Nser/Npar; %check this one
%Rpvsh_total2=Rpv_sh*Nser2/Npar2;
