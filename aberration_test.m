% Thorlabs 18917-S02 Hastings Achromatic Triplet
close all;clear;clc;
%% Aperature and Field setup
EPD=8.0;
FOV=3.5/180*pi;
NP=129;
NH=3;
py=linspace(-EPD/2,EPD/2,NP);
hy=FOV;
px=py;
hx=0;
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

%% Primary wavefront aberration terms
[PX2D,PY2D]=meshgrid(px/(EPD/2),py/(EPD/2));
[HX2D,HY2D]=meshgrid(hx/FOV,hy/FOV);
pupil_abs=double(sqrt(PX2D.^2+PY2D.^2)<=1);
r=sqrt(PX2D.^2+PY2D.^2);
phi=atan2(PY2D,PX2D);
h=sqrt(HX2D.^2+HY2D.^2);

% w040=4.326e-3;
% W=w040*r.^4;
w131=0.01/2;
W=w131*r.^3*h.*cos(phi);
% w222=0.005;
% W=w222*r.^2*h^2.*cos(phi).^2;
% w220=1.65259e-3;
% W=w220*r.^2*h^2;
% w311=9.129e-5;
% W=w311*r.*cos(phi)*h^3;
% w111=9.129e-5;
% W=w111*r.*cos(phi)*h;
figure;
imagesc(px/(EPD/2),py/(EPD/2),W.*pupil_abs);
axis image;colormap jet;

figure;
surf(px/(EPD/2),py/(EPD/2),W.*pupil_abs);

% PSF and MTF
E_XP=pupil_abs.*exp(-1i*2*pi/(wavelength_list(ind_wavelength)*1e-3)*W);
n_pad=(1024-128)/2; % see Zemax manual FFT-PSF, use pupil sampling=128
E_XP_pad=padarray(E_XP,[n_pad,n_pad]);
PSF=fftshift(fft2(ifftshift(E_XP_pad)))/numel(E_XP);
I_PSF=abs(PSF).^2;
I_PSF=I_PSF/max(I_PSF(:));
dfx=1/(1024*2/128);
fx=(-floor(size(E_XP_pad,1)/2):1:floor(size(E_XP_pad,1)/2))*dfx;
fy=fx;
f_number=5.0;
x_img=fx*2*f_number*wavelength_list(ind_wavelength);
y_img=fy*2*f_number*wavelength_list(ind_wavelength);
figure(7);
imagesc(x_img,y_img,I_PSF);
axis image;colormap jet;

MTF=fftshift(fft2(ifftshift(I_PSF)));
MTF=MTF./max(abs(MTF(:)));
center_ind=ceil(size(MTF,1)/2);
dspatial_freq=1/(128*f_number*wavelength_list(ind_wavelength)*1e-3);
spatial_freq=(-floor(size(MTF,1)/2):1:floor(size(MTF,1)/2))*dspatial_freq;
figure(8);
plot(spatial_freq(center_ind:end),abs(MTF(center_ind:end,center_ind)));hold on;
plot(spatial_freq(center_ind:end),abs(MTF(center_ind,center_ind:end)));hold off;
xlim([0,340]);

%% Ray aberrations
hx_from_OPD=-2*f_number*diff(W,1,2)/((px(2)-px(1))/(EPD/2));
hy_from_OPD=-2*f_number*diff(W,1,1)/((py(2)-py(1))/(EPD/2));

figure;
scatter(hx_from_OPD(:),hy_from_OPD(:),2);
xlim([-0.5,0.5]);ylim([-0.5,0.5]);