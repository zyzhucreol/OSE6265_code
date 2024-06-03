% Thorlabs LBF254-040 Best Form Lens - N-BK7
close all;clear;clc;
%% Aperature and Field setup
EPD=22.86;
FOV=5/180*pi;
NP=41;
NH=3;
% py=linspace(-EPD/2,EPD/2,NP);
% hy=linspace(-FOV,FOV,NH);
% px=py;
% hx=hy;
py=0;
hy=FOV;
px=0;
hx=FOV;
[PX,PY,HX,HY]=ndgrid(px,py,hx,hy);

%% Wavelength and materials setup
wavelength_list=[0.486,0.5876,0.656];
ind_wavelength=2;
B1_bk7=1.03961212;
B2_bk7=0.231792344;
B3_bk7=1.01046945;
C1_bk7=6.00069867e-3;
C2_bk7=2.00179144e-2;
C3_bk7=1.03560653e2;
n_air=1.0;
n_bk7_list=sqrt(1+B1_bk7*wavelength_list.^2./(wavelength_list.^2-C1_bk7)...
    +B2_bk7*wavelength_list.^2./(wavelength_list.^2-C2_bk7)...
    +B3_bk7*wavelength_list.^2./(wavelength_list.^2-C3_bk7));
n_bk7=n_bk7_list(ind_wavelength);

%% Prescriptions
R1=inf;
d1=0;
R2=24.02;
d2=6.5; % N-BK7
R3=-134.6;
d3=3.631074732929857E+001;
Rimg=inf;

%% Initialize paraxial rays
h1_paraxial=[PX(:),PY(:),zeros(size(PX(:)))];
t1_paraxial=zeros(size(h1_paraxial,1),1);
u1_paraxial=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];

%% Trace paraxial rays through the system
[h2_paraxial,~,t2_paraxial]=transfer_paraxial(h1_paraxial,u1_paraxial,R1,R2,d1);
[~,u2_paraxial]=refraction_paraxial(h2_paraxial,u1_paraxial,n_air,n_bk7,R2);
[h3_paraxial,~,t3_paraxial]=transfer_paraxial(h2_paraxial,u2_paraxial,R2,R3,d2);
[~,u3_paraxial]=refraction_paraxial(h3_paraxial,u2_paraxial,n_bk7,n_air,R3);
[h4_paraxial,~,t4_paraxial]=transfer_paraxial(h3_paraxial,u3_paraxial,R3,Rimg,d3);

%% Initialize real rays
h1=[PX(:),PY(:),zeros(size(PX(:)))];
r1=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];
t1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
r1=r1./t1;

%% Trace real rays through the system
[h2,~,t2]=transfer(h1,r1,R1,R2,d1);
[~,r2]=refraction(h2,r1,n_air,n_bk7,R2);
[h3,~,t3]=transfer(h2,r2,R2,R3,d2);
[~,r3]=refraction(h3,r2,n_bk7,n_air,R3);
[h4,~,t4]=transfer(h3,r3,R3,Rimg,d3);