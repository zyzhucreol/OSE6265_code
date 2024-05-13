function [ h2,r2 ] = refraction( h1,r1,n1,n2,R )
%Refraction Exact ray-tracing of refraction from single surface
%   h1, r1: incident ray coordinates and directional vector
%   n1, n2: refraction index before and after the surface
%   R: radius of curvature of the surface
%   h2, r2: refracted ray coordinates and directional vector

x1=h1(:,1);
y1=h1(:,2);
if not(isinf(R))
    z1=R-sign(R)*sqrt(R.^2-x1.^2-y1.^2);
else
    z1=zeros(size(x1));
end

%% Coordinates of emerging beam
h2=h1;

%% Directional vector of refracted ray
% normalize input directional vector
r1=r1./sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);

grad_g=-sign(R)*[2*x1,2*y1,2*(z1-R)];
grad_g=grad_g./sqrt(grad_g(:,1).^2+grad_g(:,2).^2+grad_g(:,3).^2);

r1_dot_g=r1(:,1).*grad_g(:,1)+r1(:,2).*grad_g(:,2)+r1(:,3).*grad_g(:,3);
r1_dot_g=repmat(r1_dot_g,1,3);
r_tan=r1-r1_dot_g.*grad_g;
abs_r_tan=r_tan(:,1).^2+r_tan(:,2).^2+r_tan(:,3).^2;
r2=r_tan+sqrt(n2^2/n1^2-abs_r_tan).*grad_g;

% normalize the directional vector of refracted ray
r2=r2./sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);

end