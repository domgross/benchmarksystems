s=tf('s');

T_s_sim=1e-7;%1e-5;
T_s_cont=1e-4;%1e-4;

f_b=50;

w_b=2*pi*f_b;

%Grid base values
omega_base=2*pi*f_b;
S_LV_base=100e6;
V_HV_base=230e3; %L-L rms
V_MV_base=30e3; 
V_LV_base=3.3e3;%L-L rms
I_LV_base=S_LV_base/V_LV_base;
Z_LV_base=V_LV_base/I_LV_base;
L_LV_base=Z_LV_base/omega_base;
C_LV_base=1/(Z_LV_base*omega_base);

%GSC/RSC base values
V_dc_base=7.92e3;
S_GSC_base=S_LV_base;

I_GSC_base=S_GSC_base/V_LV_base;

%GSC RLC filter parameters
R_f_pu=0.01;
L_f_pu=0.1;
C_f_pu=0.05;

R_f=R_f_pu*Z_LV_base;
L_f=L_f_pu*L_LV_base;
C_f=C_f_pu*C_LV_base;

%DC-link capacitor
H_DC_base=150e-3;
C_dc=2*S_LV_base*H_DC_base/V_dc_base^2;

%WT and PMSG 
S_WTPMSG_base=20*5e6; %rated power of a single WT
J_PMSGWT_LS=38759228; %Kg m^2 at rotor speed of 12.1 rpm
D_PMSGWT_LS=6215e3; %Nm/(rad/s) at rotor speed of 12.1 rpm

RotorRPM=12.1;
MachineRPM=1173.7;
H_PMSGWT_LS=J_PMSGWT_LS*(RotorRPM/60*2*pi)^2/(2*S_WTPMSG_base);

%usable energy
Ediff=H_PMSGWT_LS-J_PMSGWT_LS*(RotorRPM/60*2*pi*0.99)^2/(2*S_WTPMSG_base);

S_PMSG_base=S_WTPMSG_base;
J_PMSG_WT_MS=2*H_PMSGWT_LS*S_PMSG_base/(2*pi*f_b)^2; %kg m^2 at electrical speed
D_PMSG_WT_MS=D_PMSGWT_LS*2*pi/(MachineRPM/RotorRPM*3)^2;

L_PMSG_PU=0.1;
R_PMSG_PU=0.05;

L_PMSG=L_PMSG_PU*L_LV_base;
R_PMSG=R_PMSG_PU*Z_LV_base;

%Initialization
RSC_theta0=180/180*pi;
GSC_theta0=-45.0572/180*pi*0-60/180*pi;
PMSG_theta0=15/180*pi;

j=[0,-1;1,0];
Z_PMSG=R_PMSG*eye(2)+L_PMSG*w_b*j;
T_park=@(theta) 2/3*[sin(theta),sin(theta-2/3*pi),sin(theta+2/3*pi);
                     cos(theta),cos(theta-2/3*pi),cos(theta+2/3*pi);
                     1/2,1/2,1/2];
    
I_PMSG_0=T_park(0)^-1*[Z_PMSG^-1*V_LV_base*sqrt(2/3)*([cos(RSC_theta0);sin(RSC_theta0)]-[cos(PMSG_theta0);sin(PMSG_theta0)]);0];

%Pitch motor and controller
T_pitch_servo=0.5; %seconds
PI_pitch.Ka=10;
PI_pitch.Tr=5;
PI_pitch.Ta=1;

%wind speed in m/s
v_wind=12;

%Turbine radius in m
r_wt=45;

%Aerodynamic time constant according to Knudsen and Bak, Simple Model for Describing and Estimating Wind Turbine Dynamic Inflow, 2013: 
T_wt=3*r_wt/v_wind;

%First order approximation of the aerodynamic time constant
G_aero=ss(c2d(1/(T_wt*s+1),T_s_sim));

%Aerodynamic model guessed from data (overshoot of factor 2, settling time
%between 5-10 seconds)
G_aero=ss(c2d(2-1/(2*s+1),T_s_sim));


G_pitch=PI_pitch.Ka*(PI_pitch.Tr*1/s+1)/(PI_pitch.Ta*s);
G_pitch_disc=ss(c2d(G_pitch,T_s_cont,'tustin'));

beta_0=0;%10; %10; for value different than 0, Kp for pitch blade controller has to be included
Kp_beta=270*0; % 0 for beta0=0;

omega_c=omega_base/10;
G_lowpass=ss(1/(1/omega_c*s+1));
G_lowpass_disc=ss(c2d(G_lowpass,T_s_cont,'tustin'));

w0_f1=2*pi*f_b/10;
w20_f1=2*pi*2*f_b;
Q_f1=3;
s=tf('s');
Notch_w_f1=ss(c2d((s^2+w0_f1^2)/(s^2+w0_f1/Q_f1*s+w0_f1^2),T_s_cont,'tustin'));
Notch_w_f1_cont=ss((s^2+w0_f1^2)/(s^2+w0_f1/Q_f1*s+w0_f1^2));

omega_lp1khz=2*pi*1e3*2;
G_lowpass_1khz=ss(1/(1/omega_lp1khz*s+1));

G_lowpass_sim_1khz=ss(c2d(G_lowpass_1khz,T_s_sim,'tustin'));

G_lowpass_disc_1khz=ss(c2d(G_lowpass_1khz,T_s_cont,'tustin'));

omega_lp_help=1*2*pi;
G_lowpass_help=ss(1/(1/omega_lp_help*s+1));
G_lowpass_sim_help=ss(c2d(G_lowpass_help,T_s_sim,'tustin'));


w0_f10khz=2*pi*10e3;
Q_f10khz=5;
Notch_w_10khz=ss(c2d((s^2+w0_f10khz^2)/(s^2+w0_f10khz/Q_f10khz*s+w0_f10khz^2),T_s_sim,'tustin'));
Notch_w_10khz_cont=ss(c2d((s^2+w0_f10khz^2)/(s^2+w0_f10khz/Q_f10khz*s+w0_f10khz^2),T_s_cont,'tustin'));


omega_wt=1.4;
omega_wt_r=1.435; %was 1.4

k_theta_GSC=0.1;
k_p_GSC=0.02; %0.001*25; %0.02; %(this was was always used )%0.015;

k_theta_MSC=2/0.6; %5; %put 6 instead of 5 %0.5 try
k_p_MSC=k_p_GSC;%*k_theta_MSC/k_theta_GSC;


G_eb_c=(1+0.002*s)/(3e-4*s+1)+20/(2*s+1);
G_eb_MSC_c=(1+0.025*s)/(1e-3*s+1);
G_eb_GSC_c=k_theta_GSC+(k_p_GSC*s)/(1e-3*s+1);
G_diff_cont=ss(s/(20e-3/3*s+1));
G_diff=ss(c2d(s/(20e-3/3*s+1),T_s_cont,'tustin'));

G_eb_MSC=ss(c2d(G_eb_MSC_c,T_s_cont,'tustin'));
G_eb_GSC=ss(c2d(G_eb_GSC_c,T_s_cont,'tustin'));
