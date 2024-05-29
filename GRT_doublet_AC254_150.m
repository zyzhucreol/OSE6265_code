% Thorlabs 18917-S02 Hastings Achromatic Triplet
close all;clear;clc;
%% Aperature and Field setup
EPD=22.9;
FOV=5/180*pi;
NP=129;
NH=3;
py=linspace(-EPD/2,EPD/2,NP);
hy=linspace(-FOV,FOV,NH);
px=py;
hx=hy;
[PX,PY,HX,HY]=ndgrid(px,py,hx,hy);

%% Wavelength and materials setup
wavelength_list=[0.4,0.7,1.1];
ind_wavelength=2;
n_air=1.0;
B1_LAK22=1.14229781;
B2_LAK22=0.535138441;
B3_LAK22=1.04088385;
C1_LAK22=5.85778594e-3;
C2_LAK22=1.98546147e-2;
C3_LAK22=1.00834017e2;
n_LAK22_list=sqrt(1+B1_LAK22*wavelength_list.^2./(wavelength_list.^2-C1_LAK22)...
    +B2_LAK22*wavelength_list.^2./(wavelength_list.^2-C2_LAK22)...
    +B3_LAK22*wavelength_list.^2./(wavelength_list.^2-C3_LAK22));
n_LAK22=n_LAK22_list(ind_wavelength);
B1_SF10=1.62153902;
B2_SF10=0.256287842;
B3_SF10=1.64447552;
C1_SF10=1.22241457e-2;
C2_SF10=5.95736775e-2;
C3_SF10=1.47468793e2;
n_SF10_list=sqrt(1+B1_SF10*wavelength_list.^2./(wavelength_list.^2-C1_SF10)...
    +B2_SF10*wavelength_list.^2./(wavelength_list.^2-C2_SF10)...
    +B3_SF10*wavelength_list.^2./(wavelength_list.^2-C3_SF10));
n_SF10=n_SF10_list(ind_wavelength);

%% Initialize rays
h1=[PX(:),PY(:),zeros(size(PX(:)))];
r1=[tan(HX(:)),tan(HY(:)),ones(size(HX(:)))];
t1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
r1=r1./t1;

%% Prescriptions
R1=inf;
d1=0;
R2=87.861;
d2=6; % N-LAK22
R3=-105.644;
d3=3; % N-SF10
R4=inf;
d4=143.899;
Rimg=inf;

%% Trace the rays through the system
[h2,~,t2]=transfer(h1,r1,R1,R2,d1);
[~,r2]=refraction(h2,r1,n_air,n_LAK22,R2);
[h3,~,t3]=transfer(h2,r2,R2,R3,d2);
[~,r3]=refraction(h3,r2,n_LAK22,n_SF10,R3);
[h4,~,t4]=transfer(h3,r3,R3,R4,d3);
[~,r4]=refraction(h4,r3,n_SF10,n_air,R4);
[h5,~,t5]=transfer(h4,r4,R4,Rimg,d4);

%% Plotting
himg=reshape(h5,length(px),length(py),length(hx),length(hy),3);

plot_ind_hx=2;
plot_ind_hy=2;

% standard spot diagram
h_ssd=reshape(himg(:,:,plot_ind_hx,plot_ind_hy,:),[],3);
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
dXP=-149.2825; % location of XP
% dXP=0; % other reference spheres located at arbitrary plane
[~,~,t_corr]=transfer(h5,r4,Rimg,-dXP,dXP);
eikonal=t1+t2+t3*n_LAK22+t4*n_SF10+t5+t_corr;
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

%% PSF and MTF
E_XP=pupil_abs.*exp(-1i*2*pi*OPD2D);
n_pad=(1024-128)/2; % see Zemax manual FFT-PSF, use pupil sampling=128
E_XP_pad=padarray(E_XP,[n_pad,n_pad]);
PSF=fftshift(fft2(ifftshift(E_XP_pad)))/numel(E_XP);
I_PSF=abs(PSF).^2;
I_PSF=I_PSF/max(I_PSF(:));
dfx=1/(1024*2/128);
fx=(-floor(size(E_XP_pad,1)/2):1:floor(size(E_XP_pad,1)/2))*dfx;
fy=fx;
f_number=6.50007;
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
xlim([0,230]);

%% Wavefront and ray aberrations
hx_from_OPD=-abs(dXP)*diff(OPD_px*(wavelength_list(ind_wavelength)*1e-3))/(px(2)-px(1));
hy_from_OPD=-abs(dXP)*diff(OPD_py*(wavelength_list(ind_wavelength)*1e-3))/(py(2)-py(1));

figure(2);
hold on;
plot(px(1:end-1),hx_from_OPD);hold off;
figure(3);
hold on;
plot(py(1:end-1),hy_from_OPD);hold off;

OPD_from_hx=-cumsum(h_rf_px(:,1)/abs(dXP))*(px(2)-px(1))/(wavelength_list(ind_wavelength)*1e-3);
OPD_from_hy=-cumsum(h_rf_py(:,2)/abs(dXP))*(py(2)-py(1))/(wavelength_list(ind_wavelength)*1e-3);
figure(4);
hold on;
plot(px,OPD_from_hx);hold off;
figure(5);
hold on;
plot(py,OPD_from_hy);hold off;