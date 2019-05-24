function [Ke]=vigorotar(Kl,n1,n2,vz)
%vz es la auxiliar. Indica direccion de z
vx=n2-n1;
vx=vx/norm(vx);

vz=vz/norm(vz);
vy =cross(vz,vx);
vy=vy/norm(vy);
lambda = [ vx ; vy ; vz ];
T=blkdiag(lambda,lambda,lambda,lambda);
Ke = T.'*Kl*T;
return
%La forma que lo hizo Battisti.
% No entiendo nada de esto

% v1 = n2-n1;
% vd1 = v1/norm(v1);
% vp = vz-n1;
% vd3 = cross(vd1,vp)/norm(cross(vd1,vp));
% vd2 = cross(vd3,vd1);
% lambda = [ vd1 ; vd2 ; vd3 ];
% T=blkdiag(lambda,lambda,lambda,lambda);
% Ke = T.'*Kl*T;
end


