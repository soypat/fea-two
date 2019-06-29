syms x
E=2;
b=1; h=2;
A=b*h;
L=.5;
Iz = b*h^3/12;
Iy = h*b^3/12;
nu=0.3;
k = 2;
G=E/(2-2*nu);
ay=12*E*Iy/(k*G*A*L^2);
az = 12*E*Iz/(k*G*A*L^2);
by=1/(1-ay);
bz = 1/(1-az);
N1 = 1-x;
N2 = x;

syms by ay bz by az 
Hv1 = by*(2*x^3-3*x^2 + ay*x+1 -ay);
Hv2 = by*(-2*x^3+3*x^2 - ay*x);
Hw1 = bz*(2*x^3 - 3*x^2+az*x+1-az);
Hw2 = bz*(-2*x^3 + 3*x^2 - az*x);
Ht1 = L*by*(x^3 + (.5*ay-2)*x^2+(1-.5*ay)*x);
Ht2 = L*by*( x^3 -(1+.5*ay)*x^2 + .5*ay*x );
Hp1 = L*bz*(x^3 +(.5*az - 2)*x^2 + (1-.5*az)*x);
Hp2 = L*bz*(  x^3-(1+.5*az)*x^2 + .5*az*x);
Gv1 = 6*by/L*(x^2 - x);
Gv2 = 6*by/L*(-x^2+x);
Gw1 = 6*bz/L*(x^2-x);
Gw2 = 6*bz/L*(-x^2+x);
Gt1 = by*(3*x^2+(ay-4)*x+1-ay);
Gt2 = by*(3*x^2 - (ay+2)*x);
Gp1 = bz*(3*x^2 + (az -4)*x +1 -az);
Gp2 = bz*(3*x^2 - (az+2)*x);
N = [N1 N2 Hv1 Hv2 Hw1 Hw2 Ht1 Ht2 Hp1 Hp2 Gv1 Gv2 Gw1 Gw2 Gt1 Gt2 Gp1 Gp2];

ddN = diff(N,x,x); % M = d^2 v/dx^2



% p = psi, ph = phi, t=tita=theta
% Para interpolar valores:

% u = N1*u1+N2*u2;
syms v1 t1 v2 t2 w1 p1 w2 p2
v = Hv1*v1 + Ht1*t1 + Hv2*v2 + Ht2*t2;
w = Hw1*w1 + Hp1*p1 + Hw2*w2 +Hp2*p2;
ddv = diff(v,x,x);
ddw = diff(w,x,x);
dddv = diff(ddv,x);
dddw = diff(ddw,x);
% ph = N1*ph1 + N2*ph2;
% t = Gv1*v1+Gt1*t1 + Gt2*t2+Gv2*v2;
% p = Gw1*w1 + Gp1*p1 + Gw2*w2 + Gp2*p2;

