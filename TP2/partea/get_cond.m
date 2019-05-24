function digitos_perdidos = get_cond(K)
% CK = max(eigvals)/min(eigvals);
N=size(K,2);
S=zeros(N);
for i=1:N
    S(i,i)=sqrt(K(i,i))^(-1);
end
Ks=S*K*S;
scaledeigval=eig(Ks);
CKscaled=max(scaledeigval)/min(scaledeigval);
digitos_perdidos = floor(log10(CKscaled));
end

