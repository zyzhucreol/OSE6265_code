function [ h2,r2,t ] = transfer( h1,r1,R1,R2,d )
%transfer Exact ray-tracing for translation from surface 1 to 2
%   h1, r1: incident ray cordinates and directional vector
%   R1, R2: radius of curvature of two surfaces
%   d: seperation between the mid-point of two surfaces
%   h2, r2: transfered ray coordinates and directional vector
%   t: Path length of the ray

x1=h1(:,1);
y1=h1(:,2);
if not(isinf(R1))
    z1=R1-sign(R1)*sqrt(R1.^2-x1.^2-y1.^2);
else
    z1=zeros(size(x1));
end

% normalize r1 to safeguide the assumption in transfer quadratic equation
r1=r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);

rx1=r1(:,1);
ry1=r1(:,2);
rz1=r1(:,3);

%% Directional vector of the transfered ray
r2=r1;

%% Coordinates of ray after transfer onto the other surface
if not(isinf(R2))
    P=x1.*rx1+y1.*ry1+(z1-d-R2).*rz1;
    Q=x1.^2+y1.^2+(z1-d-R2).^2-R2.^2;
    t=-sign(R2)*sqrt(P.^2-Q)-P;
else
    t=(d-z1)./rz1;
end

h2=h1+t.*r1-repmat([0,0,d],size(h1,1),1);

end