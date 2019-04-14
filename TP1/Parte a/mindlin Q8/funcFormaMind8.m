clear dN N dNaux X dNxx dNyy dNxy 
syms ksi eta real

X = [1 eta ksi ksi.^2 ksi.*eta eta.^2 (ksi.^2).*eta ksi.*(eta.^2)];%Cambia de Q4 a Q8
% X = [1 ksi eta ksi.*eta];
% w1 t1 t2 
Xdx = diff(X,ksi);
Xdy = diff(X,eta);

uNod = [-1 -1
        1 -1
        1 1
        -1 1
        0 -1
        1 0
        0 1
        -1 0];
Nnodporelem = length(uNod);
Ndofpornod = 3; %Placas Mindlin
Ndofporelem = Nnodporelem*Ndofpornod;
if Nnodporelem ~=length(X)
    error('Funcionalidad deberia tener longitud igual a numero de nodos por elemento para Mindlin.')
end
A = zeros(Nnodporelem,length(X));

for i=1:Nnodporelem
    ksi=uNod(i,1); eta = uNod(i,2);
    A(i,:) = double(subs(X));
end

syms ksi eta real
shapefuns = X*inv(A);
N = shapefuns;
dNx = diff(shapefuns,ksi);
dNy = diff(shapefuns,eta);

dN = [dNx;dNy];
