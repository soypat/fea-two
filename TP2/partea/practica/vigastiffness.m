function [Ke] = vigastiffness(E,nu,A,Iz,Iy,K,L)
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
end

