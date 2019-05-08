function [Ke]=vigorotar(Kl,p1,p2,p3)
% p1=[0 0 0];
% p2=[1 1 0];
% p3=[0 0 1];
v1 = p2-p1;
vd1 = v1/norm(v1);
vp = p3-p1;
vd3 = cross(vd1,vp)/norm(cross(vd1,vp));
vd2 = cross(vd3,vd1);
lambda = [ vd1 ; vd2 ; vd3 ];
T=blkdiag(lambda,lambda,lambda,lambda);
Ke = T.'*Kl*T;
end


