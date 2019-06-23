function [T] = Tv(vx,sz)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
vy =cross(sz,vx);
vz =cross(vx,vy);
lambda = [ vx/norm(vx) ; vy/norm(vy) ; vz/norm(vz) ];
T=blkdiag(lambda,lambda,lambda,lambda);
end

