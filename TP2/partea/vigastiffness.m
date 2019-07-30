function [Ke,Me] = vigastiffness(E,nu,rho,A,Iz,Iy,K,L)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    X=A*E/L;
    G = E/(2+2*nu);
    Y4=2*E*Iz/L;
    Y3=Y4*2;
    Y2=3*Y4/L;
    Y1=2*Y2/L;
    Z4=2*E*Iy/L;
    Z3=Z4*2;
    Z2=3*Z4/L;
    Z1=2*Z2/L;
    S=G*K/L;
    Ke=[X 0 0 0 0 0 -X 0 0 0 0 0
        0 Y1 0 0 0 Y2 0 -Y1 0 0 0 Y2
        0 0 Z1 0 -Z2 0 0 0 -Z1 0 -Z2 0
        0 0 0 S 0 0 0 0 0 -S 0 0
        0 0 -Z2 0 Z3 0 0 0 Z2 0 Z4 0
        0 Y2 0 0 0 Y3 0 -Y2 0 0 0 Y4
        -X 0 0 0 0 0 X 0 0 0 0 0
        0 -Y1 0 0 0 -Y2 0 Y1 0 0 0 -Y2
        0 0 -Z1 0 Z2 0 0 0 Z1 0 Z2 0
        0 0 0 -S 0 0 0 0 0 S 0 0
        0 0 -Z2 0 Z4 0 0 0 Z2 0 Z3 0
        0 Y2 0 0 0 Y4 0 -Y2 0 0 0 Y3];
    D=0;
    a=L;
    rx=sqrt(Iz/A);
    Mediag = blkdiag(70,78,78,70*rx^2,8*a^2, 8*a^2, 70, 78, 78, 70*rx^2,8*a^2,8*a^2);
    Me =  [D 0 0 0 0 0 35 0 0 0 0 0;
                        0 D 0 0 0 22*a 0 27 0 0 0 -13*a;
                        0 0 D 0 -22*a 0 0 0 27 0 13*a 0;
                        0 0 0 D 0 0 0 0 0 -35*rx^2 0 0;
                        0 0 0 0 D 0 0 0 -13*a 0 -6*a^2 0;
                        0 0 0 0 0 D 0 13*a 0 0 0 -6*a^2;
                        0 0 0 0 0 0 D 0 0 0 0 0;
                        0 0 0 0 0 0 0 D 0 0 0 -22*a;
                        0 0 0 0 0 0 0 0 D 0 22*a 0;
                        0 0 0 0 0 0 0 0 0 D 0 0;
                        0 0 0 0 0 0 0 0 0 0 D 0;
                        0 0 0 0 0 0 0 0 0 0 0 D];
   Me = rho*A*L/210 *(Me+Me'+Mediag);
end

