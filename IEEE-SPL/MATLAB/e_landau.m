function thetas = e_landau(z)
%E_LANDAU Computes the Exact Landau solution for circle fitting
%   THETAS = e_landau(Z) where Z is a N-by-2 matrix of N 2D vectors 
%   distributed on a circle returns a vector of size 3 where the first 
%   component is the estimated radius and the two others are the estimated
%   center.
%
%   The algorithm is given in S. Thomas and A. Chan, “Simple Approach for 
%   the Estimation of Circular Arc Center and Its Radius,” Computer Vision,
%   Graphics, and   Image Processing, vol. 45, pp. 362–370, Mar. 1989. 

n=size(z,1);

x = z(:,1);
y = z(:,2);

a1 = 2*(sum(x)^2-n*sum(x.^2));
b1 = 2*(sum(x)*sum(y)-n*sum(x.*y));
b2 = 2*(sum(y)^2-n*sum(y.^2));
c1 = (sum(x.^2)*sum(x)-n*sum(x.^3)+sum(x)*sum(y.^2)-n*sum(x.*y.^2));
c2 = (sum(x.^2)*sum(y)-n*sum(y.^3)+sum(y)*sum(y.^2)-n*sum(x.^2.*y));

xc = (c1*b2-c2*b1)/(a1*b2-b1^2);
yc = (a1*c2-b1*c1)/(a1*b2-b1^2);
r = sqrt(1/n*(sum(x.^2)-2*xc*sum(x)+n*xc^2+sum(y.^2)-2*yc*sum(y)+n*yc^2));

thetas = [r; xc; yc];

end

