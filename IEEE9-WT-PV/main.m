close all; clear; clc;

%% Wind turbine parameters based on aggregate of 5MW wind turbines
data_WT

%% PV parameters for ATU Optronics PM060MBR_305W from Simscape
data_PV

%% HVDC parrameters
data_HVDC

%% Network base values
S_b=100*(10^6);
V_b=230*(10^3);
f_b=50;
w_b=2*pi*f_b;
P_b=S_b;Q_b=S_b;
I_b=S_b/(sqrt(3)*V_b);
Z_b=(V_b^2)/S_b;
L_b=Z_b/w_b;
C_b=1/(w_b*Z_b);

%% Branch parameters 
R_lines=[0,0.017,0.039,0,0.0119,0.0085,0,0.032,0.01]*Z_b;
L_lines=[0.0576,0.092,0.17,0.0586,0.1008,0.072,0.0625,0.161,0.085]*L_b;
R_14=R_lines(1); L_14=L_lines(1);
R_45=R_lines(2); L_45=L_lines(2);
R_56=R_lines(3); L_56=L_lines(3);
R_36=R_lines(4); L_36=L_lines(4);
R_67=R_lines(5); L_67=L_lines(5);
R_78=R_lines(6); L_78=L_lines(6);
R_82=R_lines(7); L_82=L_lines(7);
R_89=R_lines(8); L_89=L_lines(8);
R_94=R_lines(9); L_94=L_lines(9);

%% Shunt parameters ([pu]*bas≠e)
C_4=0.1670*C_b;C_6=0.2835*C_b;C_8=0.2275*C_b;

%% Converter params
V1_rms=1000; %Low voltage side l-l
I_b_LV=S_b/(sqrt(3)*V1_rms);
V_m=sqrt(2/3)*V1_rms;
Vdc_n=3*V_m;
n=200;
C_dc=0.008*(n);


%Converter base values
V_C_base=V1_rms;
omega_base=2*pi*50;
S_C_base=S_b/3;
I_C_base=S_C_base/V_C_base;
Z_C_base=V_C_base/I_C_base;
L_C_base=Z_C_base/omega_base;
C_C_base=1/(Z_C_base*omega_base);
L_f_pu=0.1;
R_f_pu=L_f_pu/10;
C_f_pu=0.05;

L_f=L_f_pu*L_C_base;
R_f=R_f_pu*Z_C_base;
C_f=C_f_pu*C_C_base;

L_g_pu=L_f_pu/3;
R_g_pu=L_g_pu/10;

L_g=L_g_pu*L_C_base;
R_g=R_g_pu*Z_C_base;

R_dc=(Vdc_n/(0.05*(S_b)/Vdc_n));
w_f=2*pi*5;

%hvdc parameters
Vhvdc_dc=320e3;
Vhvdc_m=Vhvdc_dc/3;% phase to groung 
Vhvdc_rms=Vhvdc_m*sqrt(3/2);  %phase to phase rms
Shvdc_b=210e6; %900e6;

Ihvdc_b_dc=Shvdc_b/Vhvdc_dc;
Ihvdc_b_LV=Shvdc_b/(3/sqrt(2)*Vhvdc_m);


%% LV/MV transformer parameters
m=100;
S_bt=S_b;
V2_rms=13800;%Medium voltage side
R1_pu=1*0.00734/m;
L1_pu=1*0.0186/m;
R2_pu=R1_pu;
L2_pu=L1_pu;
Rm_pu=10*347.82/m; %1^2/13.8^2*100/1.6/100*347.82; %10*347.82/m;
Lm_pu=10*156.68/m;%1^2/13.8^2*100/1.6/100*156.68;%10*156.68/m;
%% PV base parametrs (needed for pu->si controller gains)
Vpv_base=Vdc_n; %Vpv_mppmax;
Vpv_m=Vpv_base/3;
Ipv_b_LV=S_b/(3/sqrt(2)*Vpv_m);

%% Control parameters

%DC source and governor-turbine time constants
tau_dc=0.05;
tau_g=5;
G_cont_sm=ss(-1/tau_g,1/tau_g,1,[]);
G_sm_dis=ss(c2d(G_cont_sm,T_s_sim,'tustin'));

%defining SM governer gain----------------------

fdr_m=0.05;

% grid-forming converter control----------------
I_b_dc=S_b/Vdc_n;


%P-f droop
fdr=0.05; 
fdr_p=0.001; 
m_p=fdr_p/Vdc_n*w_b; 
mpv_p=fdr_p/Vdc_n*w_b; 
mpv_p1=fdr_p/Vdc_n*w_b;
mpv_p2=fdr_p/Vdc_n*w_b;
mpv_p3=fdr_p/Vdc_n*w_b;
mpv_p4=fdr_p/Vdc_n*w_b;
mhvdc_p=0.001/Vhvdc_dc*w_b; 

mhvdc_p1=0.001/Vhvdc_dc*w_b;
mhvdc_p2=0.001/Vhvdc_dc*w_b;

%Q-V droop
fdr_q=0.01;
m_q=fdr_q*V_m/S_b;
mpv_q=fdr_q*Vpv_m/S_b;
mpv_q1=fdr_q*Vpv_m/S_b;
mpv_q2=fdr_q*Vpv_m/S_b;
mpv_q3=fdr_q*Vpv_m/S_b;
mpv_q4=fdr_q*Vpv_m/S_b;

mhvdc_q=0.01*Vhvdc_m/Shvdc_b;

mhvdc_q1=0.01*Vhvdc_m/Shvdc_b;

mhvdc_q2=0.01*Vhvdc_m/Shvdc_b;

%Vdc-f droop
fdc=1/5;
m_dc=fdr*4*w_b/Vdc_n;


mpv_dc=fdc*w_b/Vdc_n;
kg1=-(1.24665-1.30193)/(1.30256-1.2801);
kg1_r=4.2408; 
kg2=1;
kg2_r=2.1133;

mpv_dc1=fdr*kg1*w_b/Vdc_n; 
mpv_dc2=fdr*kg2*w_b/Vdc_n;

mpv_dc1_r=fdr*kg1_r*w_b/Vdc_n;
mpv_dc2_r=fdr*kg2_r*w_b/Vdc_n;

kgwt=0.6;
Pwt_max=0.75*S_b;

mpv_dc3=fdc*w_b/Vdc_n;
mpv_dc4=fdc*w_b/Vdc_n;

mhvdc_dc=fdc*w_b/Vhvdc_dc; %was fdc
mhvdc_dc2=fdc*w_b/Vhvdc_dc*1.75;%*1.75;

kg3_r=2.6780;
kg4_r=2.7517;

mpv_dc3_r=fdr*kg3_r*w_b/Vdc_n;
mpv_dc4_r=fdr*kg4_r*w_b/Vdc_n;

%Vdc - Idc droop
k_dc=S_b/Vdc_n*fdc/fdr*1/Vdc_n;

s=tf('s');
G=w_f/(s+w_f);
pfc=ss(G);
pfd=ss(c2d(G,T_s_cont));
pf0=0;


s=tf('s');
G=1/(.1*s+1);
dfc=ss(G);
df0=0;     

%high pass virtual resistor
T_hp=1e-1;
R_virt=0.4;
G_hp=ss(R_virt*eye(3)*(1-1/(T_hp*s+1)));
G_hp_disc=ss(c2d(G_hp,T_s_cont));

%terminal voltage PI
kpvpi=0.4;
kivpi=1;
kinitpi=1;
kffpi=0;

%HVDC voltage PI
kpvpi_hvdc=4.5;
kivpi_hvdc=10;
kffpi_hvdc=0;

kpvpi_WT=0.4;
kivpi_WT=1;
kinitpi_WT=1;

%k_dc=eta_1/(Vdc_n*m_p);
K_p=(1/Vdc_n);%(1/Vdc_n);
K_r=1/R_dc;

Kp_v =0.15; %0.15(good) %0.07;%2.3; %1.9 (dominic)%7.3;%7.57; %7.5*I_b_LV/V_m; %in pu
Ki_v =0.69;%0.69(good) %0.1; %0.83 %0.69(dominic) %0.75; %0.97; %1.13*I_b_LV/V_m; %in pu
Kff_v = 1;
       
Kp_i =2.1 %2.1 %1.7 %2.1; %2.7 looked very nice
Ki_i = 0.79;
Kff_i = 1;

K_WT_p_v =0.15;%0.15 %0.07;%2.3; %1.9 (dominic)%7.3;%7.57; %7.5*I_b_LV/V_m; %in pu
K_WT_i_v =0.69;%0.69 %0.1; %0.83 %0.69(dominic) %0.75; %0.97; %1.13*I_b_LV/V_m; %in pu
K_WT_ff_v = 1;
       
K_WT_p_i =2.1  %1.3 %1.3
K_WT_i_i = 0.79;
K_WT_ff_i = 1;


%%%%     
     % Voltage loop----------------------------------------
     Kvscp_v = Kp_v*I_b_LV/V_m; %7.5*I_b_LV/V_m;  %0.52
     Kvsci_v = Ki_v*I_b_LV/V_m; %1.13*I_b_LV/V_m; 
     Kff_v = 1;
        
     % Current loop
     Kvscp_i = Kp_i*V_m/I_b_LV;%0.074*V_m/I_b_LV; %0.04*V_m/I_b_LV;
     Kvsci_i = Ki_i*V_m/I_b_LV; 
     Kff_i = 1;


     %Voltage loops PV
Kpvp_v = Kp_v*Ipv_b_LV/Vpv_m; %2.91*Ipv_b_LV/Vpv_m; %7.5*Ipv_b_LV/Vpv_m;
Kpvi_v = Ki_v*Ipv_b_LV/Vpv_m; %0.97*Ipv_b_LV/Vpv_m; %1.13
Kffpv_v = 1;
         
 % Current loop PV
Kpvp_i = Kp_i*Vpv_m/Ipv_b_LV; %0.27*Vpv_m/Ipv_b_LV; %0.04*Vpv_m/Ipv_b_LV;
Kpvi_i = Ki_i*Vpv_m/Ipv_b_LV; %1.03*Vpv_m/xIpv_b_LV;
Kffpv_i = 1;


%Voltage loops PV1
Kpvp_v1 = Kp_v*Ipv_b_LV/Vpv_m; %2.91*Ipv_b_LV/Vpv_m; %7.5*Ipv_b_LV/Vpv_m;
Kpvi_v1 = Ki_v*Ipv_b_LV/Vpv_m; %0.97*Ipv_b_LV/Vpv_m; %1.13
        
 % Current loop PV1
Kpvp_i1 = Kp_i*Vpv_m/Ipv_b_LV; %0.27*Vpv_m/Ipv_b_LV; %0.04*Vpv_m/Ipv_b_LV;
Kpvi_i1 = Ki_i*Vpv_m/Ipv_b_LV; %1.03*Vpv_m/Ipv_b_LV;

%Voltage loops PV2
Kpvp_v2 = Kp_v*Ipv_b_LV/Vpv_m; %2.91*Ipv_b_LV/Vpv_m; %7.5*Ipv_b_LV/Vpv_m;
Kpvi_v2 = Ki_v*Ipv_b_LV/Vpv_m; %0.97*Ipv_b_LV/Vpv_m; %1.13
        
 % Current loop PV2
Kpvp_i2 = Kp_i*Vpv_m/Ipv_b_LV; %0.27*Vpv_m/Ipv_b_LV; %0.04*Vpv_m/Ipv_b_LV;
Kpvi_i2 = Ki_i*Vpv_m/Ipv_b_LV; %1.03*Vpv_m/Ipv_b_LV;


s=tf('s');
GvscF=1/(L_f*s+0*R_f);
Gmod=1/(1/(2*10e3)*s+1);
GvscC=1/(C_f*s);
GvscPI_i=Kvscp_i+Kvsci_i/s;
GvscPI_v=Kvscp_v+Kvsci_v/s;

Gcl_vsc_i=(GvscPI_i*Gmod*GvscF)/(1+GvscPI_i*Gmod*GvscF);
Gcl_vsc_v=(GvscPI_v*GvscC)/(1+GvscPI_v*GvscC);
Gcl_vsc_iv=(GvscPI_v*Gcl_vsc_i*GvscC)/(1+GvscPI_v*Gcl_vsc_i*GvscC);

figure(1),subplot(2,1,1),hold on
step(Gcl_vsc_i,0.006)
step(Gcl_vsc_v,0.006)
grid on
legend('closed-loop current dynamics','approx. closed-loop voltage dynamics')
subplot(2,1,2)
step(Gcl_vsc_iv,0.06)
grid on
legend('closed-loop voltage dynamics')

%Vdc - Idc droop

%voltage loop hvdc
Khvdcp_v = 0.15*Ihvdc_b_LV/Vhvdc_m; 
Khvdci_v = 0.69*Ihvdc_b_LV/Vhvdc_m; 
Kff_hvdc_i = 1;
        
      
% Current loop
Khvdcp_i = 2.1*Vhvdc_m/Ihvdc_b_LV; 
Khvdci_i = 0.79*Vhvdc_m/Ihvdc_b_LV; 
Kff_hvdc_v = 1;
T_HVDC_PI_en=0;

%voltage loop hvdc
Khvdcp_v1 = 0.15*Ihvdc_b_LV/Vhvdc_m; 
Khvdci_v1 = 0.69*Ihvdc_b_LV/Vhvdc_m;
      
% Current loop
Khvdcp_i1 = 2.1*Vhvdc_m/Ihvdc_b_LV; 
Khvdci_i1 = 0.79*Vhvdc_m/Ihvdc_b_LV; 
s=tf('s');
GhvdcF=1/(Lhvdc_f*s+0*Rhvdc_f);
GhvdcC=1/(Chvdc_f*s);
GhvdcPI_i=Khvdcp_i+Khvdci_i/s;
GhvdcPI_v=Khvdcp_v+Khvdci_v/s;

Gcl_i=(GhvdcPI_i*Gmod*GhvdcF)/(1+GhvdcPI_i*Gmod*GhvdcF);
Gcl_v=(GhvdcPI_v*GhvdcC)/(1+GhvdcPI_v*GhvdcC);
Gcl_iv=(GhvdcPI_v*Gcl_i*GhvdcC)/(1+GhvdcPI_v*Gcl_i*GhvdcC);

figure(2),subplot(2,1,1),hold on
step(Gcl_i,0.006)
step(Gcl_v,0.006)
grid on
legend('closed-loop current dynamics','approx. closed-loop voltage dynamics')
subplot(2,1,2)
step(Gcl_iv,0.006)
grid on
legend('closed-loop voltage dynamics')

GhvdcPI_i1=Khvdcp_i1+Khvdci_i1/s;
GhvdcPI_v1=Khvdcp_v1+Khvdci_v1/s;

Gcl_i1=(GhvdcPI_i1*Gmod*GhvdcF)/(1+GhvdcPI_i1*Gmod*GhvdcF);
Gcl_v1=(GhvdcPI_v1*GhvdcC)/(1+GhvdcPI_v1*GhvdcC);
Gcl_iv1=(GhvdcPI_v1*Gcl_i*GhvdcC)/(1+GhvdcPI_v1*Gcl_i*GhvdcC);



%% Load and IBR set-points
base=2.25; % base load
load_step=0.075*S_b;% load disturbance
pl_main=S_b*base/3; %loads in [W]

main_p= 0.575;
main_p_r=1.2;

island_p=[0.68741 1.2838];
island_p_r=[0.5791 1.1923];

island_Vdc=[1 1.2853    Vpv_mppmax2/Vdc_n];
island_Vdc_r=[1 1.3430 1.4716];
island_Vm=[1 1 1];
Q_gsc_ref=0;
tie_power=0.5693;

%%

island_q=[-0.1312   -0.0230   -0.0343];
island_q_r=[-0.1409   -0.0622   -0.0525];



%set_points
%%
tau=1e-5;
%%
%
T_ms=0.001; %sampling time for output data
%%
G_setpoint_lowpass_disc=ss(c2d(ss(1/(0.1*s+1)),T_s_cont,'tustin'));
G_setpoint_lowpass_sim=ss(c2d(ss(1/(0.1*s+1)),T_s_sim,'tustin'));

B1=40*10^6;
B2=10*10^6;
TN1=50;
Cp1=0.03;
TN2=TN1;
Cp2=Cp1;
T_second=50;
pl_island=pl_main;

Tend=50;
T_redispatch=25;
T_load=20;
count=00;

load_step1=load_step;
load_step2=load_step;
data=strcat('data',num2str(count),'.mat');
