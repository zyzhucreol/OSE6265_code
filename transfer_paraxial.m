function [ h2,r2,t ] = transfer_paraxial( h1,r1,R1,R2,d )
%transfer Exact ray-tracing for translation from surface 1 to 2
%   h1, r1: incident ray cordinates and directional vector
%   R1, R2: radius of curvature of two surfaces
%   d: seperation between the mid-point of two surfaces
%   h2, r2: transfered ray coordinates and directional vector
%   t: Path length of the ray

%% Directionof the transfered ray
r2=r1;

%% Coordinates of ray after transfer onto the other surface
h2=h1+d*r1-repmat([0,0,d],size(h1,1),1);

t=sqrt((h2(:,1)-h1(:,1)).^2+(h2(:,2)-h1(:,2)).^2+d^2);

end