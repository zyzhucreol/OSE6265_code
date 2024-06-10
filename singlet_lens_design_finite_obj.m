% Thorlabs LBF254-040 Best Form Lens - N-BK7
close all;clear;clc;
%% Aperature and Field setup
EPD=22.86;
FOV=5/180*pi;
NP=41;
NH=3;
py=linspace(-EPD/2,EPD/2,NP);
hy=linspace(-FOV,FOV,NH);
px=py;
hx=hy;
% py=0;
% hy=FOV;
% px=0;
% hx=FOV;
[PX,PY,HX,HY]=ndgrid(px,py,hx,hy);
ind_marginal=sub2ind(size(PX),21,41,2,2);
ind_chief=sub2ind(size(PX),21,21,2,3);

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
dobj=130;
R1=inf;
d1=0;
R2=6.001999999999995E+001;
d2=4.0; % N-BK7
R3=-3.533000000000053E+002;
d3=9.733998017553688E+001;
Rimg=inf;

%% Constraints on R1 and R2
K=1/100;
c1=linspace(0,0.04,128);
c2=(K/(n_bk7-1)-c1)./(c1*d2*(n_bk7-1)/n_bk7-1);

S1=zeros(size(c1));
S2=zeros(size(c1));
S3=zeros(size(c1));
for ii=1:length(c1)
    R2=1/c1(ii);
    R3=1/c2(ii);
%% Initialize paraxial rays
h1_paraxial=[PX(:),PY(:),zeros(size(PX(:)))];
t1_paraxial=zeros(size(h1_paraxial,1),1);
% u1_paraxial=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];
Pobj=[-dobj*tan(HX(:)),-dobj*tan(HY(:)),-dobj*ones(size(HX(:)))];
r1=[PX(:)-Pobj(:,1),PY(:)-Pobj(:,2),-Pobj(:,3)];
u1_paraxial=r1./repmat(-Pobj(:,3),[1,3]);

%% Trace paraxial rays through the system
[h2_paraxial,~,t2_paraxial]=transfer_paraxial(h1_paraxial,u1_paraxial,R1,R2,d1);
[~,u2_paraxial]=refraction_paraxial(h2_paraxial,u1_paraxial,n_air,n_bk7,R2);
[h3_paraxial,~,t3_paraxial]=transfer_paraxial(h2_paraxial,u2_paraxial,R2,R3,d2);
[~,u3_paraxial]=refraction_paraxial(h3_paraxial,u2_paraxial,n_bk7,n_air,R3);
[h4_paraxial,~,t4_paraxial]=transfer_paraxial(h3_paraxial,u3_paraxial,R3,Rimg,d3);

%% Calculate primary wavefront aberration coefficients
h1m=h2_paraxial(ind_marginal,2);
h1c=h2_paraxial(ind_chief,2);
u1m=u1_paraxial(ind_marginal,2);
u1c=u1_paraxial(ind_chief,2);
A1m=n_air*(h1m/R2+u1m);
A1c=n_air*(h1c/R2+u1c);
h2m=h3_paraxial(ind_marginal,2);
h2c=h3_paraxial(ind_chief,2);
u2m=u2_paraxial(ind_marginal,2);
u2c=u2_paraxial(ind_chief,2);
A2m=n_bk7*(h2m/R3+u2m);
A2c=n_bk7*(h2c/R3+u2c);
H1=h1m*A1c-h1c*A1m;
H2=h2m*A2c-h2c*A2m;

u3m=u3_paraxial(ind_marginal,2);
u3c=u3_paraxial(ind_chief,2);

S1(ii)=-A1m^2*h1m*(u2m/n_bk7-u1m/n_air)-A2m^2*h2m*(u3m/n_air-u2m/n_bk7);
S2(ii)=-A1m*A1c*h1m*(u2m/n_bk7-u1m/n_air)-A2m*A2c*h2m*(u3m/n_air-u2m/n_bk7);
S3(ii)=-A1c^2*h1m*(u2m/n_bk7-u1m/n_air)-A2c^2*h2m*(u3m/n_air-u2m/n_bk7);

end

figure;
plot(c1,c2);

figure;
plot(c1,S1);hold on;
plot(c1,S2);hold on;
plot(c1,S3);hold off;
legend('S1','S2','S3');