function [outputArg1,outputArg2] = vigaVM(dlocal,samples)

x = 0:1/samples:1;
Nsamples = length(x);

Nx = A*E/L*(u2-u1); 
T = G*K*(tx2-tx1)/L;
My = E*Iy*ddw;
Mz = E*Iz*ddv;
Vy = E*Iz*dddv;
Vz = E*Iy*dddw;

ay=12*E*Iy/(k*G*A*L^2);
az = 12*E*Iz/(k*G*A*L^2);
by=1/(1-ay);
bz = 1/(1-az);

ddv = by*v1*(12*x - 6) - by*v2*(12*x - 6) - (by*t2*(ay - 6*x + 2))/2 + (by*t1*(ay + 6*x - 4))/2;
ddw = bz*w1*(12*x - 6) - bz*w2*(12*x - 6) - (bz*p2*(az - 6*x + 2))/2 + (bz*p1*(az + 6*x - 4))/2;
dddv = 3*by*t1 + 3*by*t2 + 12*by*v1 - 12*by*v2;
dddw = 3*bz*p1 + 3*bz*p2 + 12*bz*w1 - 12*bz*w2;
u = zeros(Nsamples,1);
v = zeros(Nsamples,1);
w = zeros(Nsamples,1);
Hv1 = by*(2*x^3-3*x^2 + ay*x+1 -ay);
Hv2 = by*(-2*x^3+3*x^2 - ay*x);
Hw1 = bz*(2*x^3 - 3*x^2+az*x+1-az);
Hw2 = bz*(-2*x^3 + 3*x^2 - az*x);
Ht1 = L*by*(x^3 + (.5*ay-2)*x^2+(1-.5*ay)*x);
Ht2 = L*by*( x^3 -(1+.5*ay)*x^2 + .5*ay*x );
Hp1 = L*bz*(x^3 +(.5*az - 2)*x^2 + (1-.5*az)*x);
Hp2 = L*bz*(  x^3-(1+.5*az)*x^2 + .5*az*x);

    


