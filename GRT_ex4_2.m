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

%% Initialize rays
h1=[PX(:),PY(:),zeros(size(PX(:)))];
r1=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];
t1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
r1=r1./t1;

%% Prescriptions
R1=inf;
d1=0;
R2=inf;
d2=5.340; % N-BK7
R3=-25.750;
d3=44.205;
Rimg=inf;

%% Trace the rays through the system
[h2,~,t2]=transfer(h1,r1,R1,R2,d1);
[~,r2]=refraction(h2,r1,n_air,n_bk7,R2);
[h3,~,t3]=transfer(h2,r2,R2,R3,d2);
[~,r3]=refraction(h3,r2,n_bk7,n_air,R3);
[h4,~,t4]=transfer(h3,r3,R3,Rimg,d3);

%% Plotting
himg=reshape(h4,length(px),length(py),length(hx),length(hy),3);

plot_ind_hx=2;
plot_ind_hy=2;

% standard spot diagram
h_ssd=reshape(himg(:,:,plot_ind_hy,plot_ind_hx,:),[],3);
figure;
scatter(h_ssd(:,1),h_ssd(:,2),8);
axis image;

% ray aberration (ray-fan)
h_rf_px=squeeze(himg(:,ceil(length(py)/2),plot_ind_hx,plot_ind_hy,:));
figure;
plot(px,h_rf_px(:,1));hold on;
plot(px,h_rf_px(:,2));hold off;
h_rf_py=squeeze(himg(ceil(length(px)/2),:,plot_ind_hx,plot_ind_hy,:));
figure;
plot(px,h_rf_py(:,1));hold on;
plot(px,h_rf_py(:,2));hold off;

% optical path difference (Reference: Absolute)
eikonal=t1+t2+t3*n_bk7+t4;
eikonal=reshape(eikonal,length(px),length(py),length(hx),length(hy));
eikonal_chief=repmat(eikonal(ceil(length(px)/2),ceil(length(py)/2),:,:),length(px),length(py),1,1);
OPD=(eikonal_chief-eikonal)/(wavelength_list(ind_wavelength)*1e-3);
OPD_px=squeeze(OPD(:,ceil(length(py)/2),plot_ind_hx,plot_ind_hy));
figure;
plot(px,OPD_px);
OPD_py=squeeze(OPD(ceil(length(px)/2),:,plot_ind_hx,plot_ind_hy));
figure;
plot(py,OPD_py);

% optical path difference (Reference: Exit Pupil)
dXP=-38.25103;
[~,~,t_corr]=transfer_XP(h4,r3,Rimg,dXP);
eikonal=t1+t2+t3*n_bk7+t4+t_corr;
eikonal=reshape(eikonal,length(px),length(py),length(hx),length(hy));
eikonal_chief=repmat(eikonal(ceil(length(px)/2),ceil(length(py)/2),:,:),length(px),length(py),1,1);
OPD=(eikonal_chief-eikonal)/(wavelength_list(ind_wavelength)*1e-3);
OPD_px=squeeze(OPD(:,ceil(length(py)/2),plot_ind_hx,plot_ind_hy));
figure;
plot(px,OPD_px);
OPD_py=squeeze(OPD(ceil(length(px)/2),:,plot_ind_hx,plot_ind_hy));
figure;
plot(py,OPD_py);