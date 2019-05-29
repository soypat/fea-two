function [Ke]=vigorotar(Kl,n1,n2,sz)
%sz es la auxiliar
vx=n2-n1;
vy =cross(sz,vx);
vz =cross(vx,vy);
lambda = [ vx/norm(vx) ; vy/norm(vy) ; vz/norm(vz) ];
T=blkdiag(lambda,lambda,lambda,lambda);
Ke = T.'*Kl*T;
return
end


