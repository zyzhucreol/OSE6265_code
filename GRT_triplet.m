% Thorlabs 18917-S02 Hastings Achromatic Triplet
close all;clear;clc;
%% Aperature and Field setup
EPD=8.0;
FOV=3.5/180*pi;
NP=129;
NH=3;
py=linspace(-EPD/2,EPD/2,NP);
hy=linspace(-FOV,FOV,NH);
px=py;
hx=hy;
[PX,PY,HX,HY]=ndgrid(px,py,hx,hy);

%% Wavelength and materials setup
wavelength_list=[0.486133,0.587562,0.656273];
ind_wavelength=2;
n_air=1.0;
B1_F2=1.39757037;
B2_F2=0.159201403;
B3_F2=1.26865430;
C1_F2=9.95906143e-3;
C2_F2=5.46931752e-2;
C3_F2=1.19248346e2;
n_F2_list=sqrt(1+B1_F2*wavelength_list.^2./(wavelength_list.^2-C1_F2)...
    +B2_F2*wavelength_list.^2./(wavelength_list.^2-C2_F2)...
    +B3_F2*wavelength_list.^2./(wavelength_list.^2-C3_F2));
n_F2=n_F2_list(ind_wavelength);
B1_K5=1.08511833;
B2_K5=0.199562005;
B3_K5=0.930511663;
C1_K5=6.61099503e-3;
C2_K5=2.41108660e-2;
C3_K5=1.11982777e2;
n_K5_list=sqrt(1+B1_K5*wavelength_list.^2./(wavelength_list.^2-C1_K5)...
    +B2_K5*wavelength_list.^2./(wavelength_list.^2-C2_K5)...
    +B3_K5*wavelength_list.^2./(wavelength_list.^2-C3_K5));
n_K5=n_K5_list(ind_wavelength);

%% Initialize rays
h1=[PX(:),PY(:),zeros(size(PX(:)))];
r1=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];
t1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
r1=r1./t1;

%% Prescriptions
R1=inf;
d1=0;
R2=30.2;
d2=2.5; % N-F2
R3=14.223;
d3=17; % N-K5
R4=-14.223;
d4=2.5; % N-F2
R5=-30.2;
d5=32.204;
Rimg=inf;

%% Trace the rays through the system
[h2,~,t2]=transfer(h1,r1,R1,R2,d1);
[~,r2]=refraction(h2,r1,n_air,n_F2,R2);
[h3,~,t3]=transfer(h2,r2,R2,R3,d2);
[~,r3]=refraction(h3,r2,n_F2,n_K5,R3);
[h4,~,t4]=transfer(h3,r3,R3,R4,d3);
[~,r4]=refraction(h4,r3,n_K5,n_F2,R4);
[h5,~,t5]=transfer(h4,r4,R4,R5,d4);
[~,r5]=refraction(h5,r4,n_F2,n_air,R5);
[h6,~,t6]=transfer(h5,r5,R5,Rimg,d5);

%% Plotting
himg=reshape(h6,length(px),length(py),length(hx),length(hy),3);

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

% optical path difference (Reference: Exit Pupil)
dXP=-50.34802; % location of XP
% dXP=0; % other reference spheres located at arbitrary plane
[~,~,t_corr]=transfer(h6,r5,Rimg,-dXP,dXP);
eikonal=t1+t2+t3*n_F2+t4*n_K5+t5*n_F2+t6+t_corr;
eikonal=reshape(eikonal,length(px),length(py),length(hx),length(hy));
eikonal_chief=repmat(eikonal(ceil(length(px)/2),ceil(length(py)/2),:,:),length(px),length(py),1,1);
OPD=(eikonal_chief-eikonal)/(wavelength_list(ind_wavelength)*1e-3);
OPD_px=squeeze(OPD(:,ceil(length(py)/2),plot_ind_hx,plot_ind_hy));
figure;
plot(px,OPD_px);
OPD_py=squeeze(OPD(ceil(length(px)/2),:,plot_ind_hx,plot_ind_hy));
figure;
plot(py,OPD_py);

OPD2D=OPD(:,:,plot_ind_hx,plot_ind_hy);
[PX2D,PY2D]=meshgrid(px,py);
pupil_abs=double(sqrt(PX2D.^2+PY2D.^2)<=EPD/2);
figure;
imagesc(px/(EPD/2),py/(EPD/2),OPD2D.*pupil_abs);
axis image;colormap jet;

% PSF and MTF
E_XP=pupil_abs.*exp(-1i*2*pi*OPD2D);
n_pad=(1024-128)/2; % see Zemax manual FFT-PSF, use pupil sampling=128
E_XP_pad=padarray(E_XP,[n_pad,n_pad]);
PSF=fftshift(fft2(ifftshift(E_XP_pad)))/numel(E_XP);
I_PSF=abs(PSF).^2;
I_PSF=I_PSF/max(I_PSF(:));
dfx=1/(1024*2/128);
fx=(-floor(size(E_XP_pad,1)/2):1:floor(size(E_XP_pad,1)/2))*dfx;
fy=fx;
f_number=5.04432;
x_img=fx*2*f_number*wavelength_list(ind_wavelength);
y_img=fy*2*f_number*wavelength_list(ind_wavelength);
figure;
imagesc(x_img,y_img,I_PSF);
axis image;colormap jet;

MTF=fftshift(fft2(ifftshift(I_PSF)));
MTF=MTF./max(abs(MTF(:)));
center_ind=ceil(size(MTF,1)/2);
dspatial_freq=1/(128*f_number*wavelength_list(ind_wavelength)*1e-3);
spatial_freq=(-floor(size(MTF,1)/2):1:floor(size(MTF,1)/2))*dspatial_freq;
figure;
plot(spatial_freq(center_ind:end),abs(MTF(center_ind:end,center_ind)));
xlim([0,340]);