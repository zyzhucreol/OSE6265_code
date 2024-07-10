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
B1_SF2=1.473431270;
B2_SF2=1.636818490E-01;
B3_SF2=1.369208990;
C1_SF2=1.090190980E-02;
C2_SF2=5.856836870E-02;
C3_SF2=1.274049330E+02;
n_SF2_list=sqrt(1+B1_SF2*wavelength_list.^2./(wavelength_list.^2-C1_SF2)...
    +B2_SF2*wavelength_list.^2./(wavelength_list.^2-C2_SF2)...
    +B3_SF2*wavelength_list.^2./(wavelength_list.^2-C3_SF2));
V_SF2=(n_SF2_list(2)-1)/(n_SF2_list(1)-n_SF2_list(3));
V_BK7=(n_BK7_list(2)-1)/(n_BK7_list(1)-n_BK7_list(3));

%% Constraints
Ktot=1/15;

% Use glass set SSK5, BAK1, FPL51
D=[1,1;...
    1/V_BK7,1/V_SF2];
K=D\[Ktot;0];
EFL1=1/K(1);
EFL2=1/K(2);