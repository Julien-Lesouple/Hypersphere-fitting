function thetas = fgfa(Z)
%FGFA Computes the fast geometric fit algorithm solution for sphere fitting
%   THETAS = fgfa(Z) where Z is a N-by-3 matrix of N 2D vectors 
%   distributed on a sphere returns a vector of size 4 where the first 
%   component is the estimated radius and the three others are the estimated
%   center.
%
%   The algorithm is given in Y. Sumith, “Fast Geometric Fit Algorithm for 
%   Sphere Using Exact Solution,” CoRR, 2015.

x = Z(:,1);
y = Z(:,2);
z = Z(:,3);

[N,dum]=size(x);
Sx = sum(x); Sy = sum(y); Sz = sum(z);
Sxx = sum(x.*x); Syy = sum(y.*y);
Szz = sum(z.*z); Sxy = sum(x.*y);
Sxz = sum(x.*z); Syz = sum(y.*z);
Sxxx = sum(x.*x.*x); Syyy = sum(y.*y.*y);
Szzz = sum(z.*z.*z); Sxyy = sum(x.*y.*y);
Sxzz = sum(x.*z.*z); Sxxy = sum(x.*x.*y);
Sxxz = sum(x.*x.*z); Syyz =sum(y.*y.*z);
Syzz = sum(y.*z.*z);
A1 = Sxx +Syy +Szz;
a = 2*Sx*Sx-2*N*Sxx;
b = 2*Sx*Sy-2*N*Sxy;
c = 2*Sx*Sz-2*N*Sxz;
d = -N*(Sxxx +Sxyy +Sxzz)+A1*Sx;
e = 2*Sx*Sy-2*N*Sxy;
f = 2*Sy*Sy-2*N*Syy;
g = 2*Sy*Sz-2*N*Syz;
h = -N*(Sxxy +Syyy +Syzz)+A1*Sy;
j = 2*Sx*Sz-2*N*Sxz;
k = 2*Sy*Sz-2*N*Syz;
l = 2*Sz*Sz-2*N*Szz;
m = -N*(Sxxz +Syyz + Szzz)+A1*Sz;
delta = a*(f*l - g*k)-e*(b*l-c*k) + j*(b*g-c*f);
xc = (d*(f*l-g*k) -h*(b*l-c*k) +m*(b*g-c*f))/delta;
yc = (a*(h*l-m*g) -e*(d*l-m*c) +j*(d*g-h*c))/delta;
zc = (a*(f*m-h*k) -e*(b*m-d*k) +j*(b*h-d*f))/delta;
R = sqrt(xc^2+yc^2+zc^2+(A1-2*(xc*Sx+yc*Sy+zc*Sz))/N);

thetas = [R;xc;yc;zc];

end

