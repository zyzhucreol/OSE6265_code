function [ h2,r2 ] = refraction_paraxial( h1,r1,n1,n2,R )
%Refraction Exact ray-tracing of refraction from single surface
%   h1, r1: incident ray coordinates and directional vector
%   n1, n2: refraction index before and after the surface
%   R: radius of curvature of the surface
%   h2, r2: refracted ray coordinates and directional vector

%% Coordinates of emerging beam
h2=h1;

%% Direction of emerging beam
r2=(n1-n2)/(n2*R)*h1+n1/n2*r1;
r2(:,3)=1;

end