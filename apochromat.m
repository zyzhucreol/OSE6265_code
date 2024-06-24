% Shape-dependent aberrations of a not-so-thin lens
close all;clear;clc;

%% Wavelength and materials setup
wavelength_list=[0.486,0.5876,0.656];
ind_wavelength=2;
B1_BK7=1.03961212;
B2_BK7=0.231792344;
B3_BK7=1.01046945;
C1_BK7=6.00069867e-3;
C2_BK7=2.00179144e-2;
C3_BK7=1.03560653e2;
n_air=1.0;
n_BK7_list=sqrt(1+B1_BK7*wavelength_list.^2./(wavelength_list.^2-C1_BK7)...
    +B2_BK7*wavelength_list.^2./(wavelength_list.^2-C2_BK7)...
    +B3_BK7*wavelength_list.^2./(wavelength_list.^2-C3_BK7));
B1_BAK1=1.123656620E+00;
B2_BAK1=3.092768480E-01;
B3_BAK1=8.815119570E-01;
C1_BAK1=6.447427520E-03;
C2_BAK1=2.222844020E-02;
C3_BAK1=1.072977510E+02;
n_BAK1_list=sqrt(1+B1_BAK1*wavelength_list.^2./(wavelength_list.^2-C1_BAK1)...
    +B2_BAK1*wavelength_list.^2./(wavelength_list.^2-C2_BAK1)...
    +B3_BAK1*wavelength_list.^2./(wavelength_list.^2-C3_BAK1));
B1_FPL51=1.029607000E+00;
B2_FPL51=1.880506000E-01;
B3_FPL51=7.364881650E-01;
C1_FPL51=5.168001550E-03;
C2_FPL51=1.666587980E-02;
C3_FPL51=1.389641290E+02;
n_FPL51_list=sqrt(1+B1_FPL51*wavelength_list.^2./(wavelength_list.^2-C1_FPL51)...
    +B2_FPL51*wavelength_list.^2./(wavelength_list.^2-C2_FPL51)...
    +B3_FPL51*wavelength_list.^2./(wavelength_list.^2-C3_FPL51));
B1_SF2=1.473431270;
B2_SF2=1.636818490E-01;
B3_SF2=1.369208990;
C1_SF2=1.090190980E-02;
C2_SF2=5.856836870E-02;
C3_SF2=1.274049330E+02;
n_SF2_list=sqrt(1+B1_SF2*wavelength_list.^2./(wavelength_list.^2-C1_SF2)...
    +B2_SF2*wavelength_list.^2./(wavelength_list.^2-C2_SF2)...
    +B3_SF2*wavelength_list.^2./(wavelength_list.^2-C3_SF2));
B1_SSK5=1.592226590E+00;
B2_SSK5=1.035207740E-01;
B3_SSK5=1.051740160E+00;
C1_SSK5=9.202846260E-03;
C2_SSK5=4.235300720E-02;
C3_SSK5=1.069273740E+02;
n_SSK5_list=sqrt(1+B1_SSK5*wavelength_list.^2./(wavelength_list.^2-C1_SSK5)...
    +B2_SSK5*wavelength_list.^2./(wavelength_list.^2-C2_SSK5)...
    +B3_SSK5*wavelength_list.^2./(wavelength_list.^2-C3_SSK5));
V_SSK5=(n_SSK5_list(2)-1)/(n_SSK5_list(1)-n_SSK5_list(3));
V_SF2=(n_SF2_list(2)-1)/(n_SF2_list(1)-n_SF2_list(3));
V_BK7=(n_BK7_list(2)-1)/(n_BK7_list(1)-n_BK7_list(3));
V_BAK1=(n_BAK1_list(2)-1)/(n_BAK1_list(1)-n_BAK1_list(3));
V_FPL51=(n_FPL51_list(2)-1)/(n_FPL51_list(1)-n_FPL51_list(3));
P_SSK5=(n_SSK5_list(1)-n_SSK5_list(2))/(n_SSK5_list(1)-n_SSK5_list(3));
P_SF2=(n_SF2_list(1)-n_SF2_list(2))/(n_SF2_list(1)-n_SF2_list(3));
P_BK7=(n_BK7_list(1)-n_BK7_list(2))/(n_BK7_list(1)-n_BK7_list(3));
P_BAK1=(n_BAK1_list(1)-n_BAK1_list(2))/(n_BAK1_list(1)-n_BAK1_list(3));
P_FPL51=(n_FPL51_list(1)-n_FPL51_list(2))/(n_FPL51_list(1)-n_FPL51_list(3));

%% Constraints
Ktot=1/400;

% Use glass set SSK5, BAK1, FPL51
% "Choice of material" matrix
D=[1,1,1;...
    1/V_SSK5,1/V_BAK1,1/V_FPL51;...
    P_SSK5/V_SSK5,P_BAK1/V_BAK1,P_FPL51/V_FPL51];
K=D\[Ktot;0;0];
EFL1=1/K(1);
EFL2=1/K(2);
EFL3=1/K(3);

% Use glass set N-SF2, N-BK7, S-FPL51 (N-PK52A)
% D=[1,1,1;...
%     1/V_BK7,1/V_SF2,1/V_FPL51;...
%     P_BK7/V_BK7,P_SF2/V_SF2,P_FPL51/V_FPL51];
% K=D\[Ktot;0;0];
% EFL1=1/K(1);
% EFL2=1/K(2);
% EFL3=1/K(3);