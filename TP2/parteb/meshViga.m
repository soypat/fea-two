function [nodos, elementos, elemDof] = meshViga(n1,n2,N,Ndofpornod,nodstart)
%MESHVIGA meshViga(n1,n2,Nelem,Ndofpornod,nodstart)
Nnodporelem = 2;
Ndofporelem = Nnodporelem*Ndofpornod;
nodos = zeros(N+1,3);
elemDof = zeros(N,Ndofpornod*Nnodporelem);
elementos = zeros(N,Nnodporelem);
nodos(1,:) = n1;
xv = linspace(n1(1),n2(1),N+1);
yv = linspace(n1(2),n2(2),N+1);
zv = linspace(n1(3),n2(3),N+1);
for e=1:N+1 %Iteracion sobre nodos
    nodos(e,:)=[xv(e),yv(e),zv(e)];
end

for e=1:N
       elementos(e,:) = [e e+1];
end
e=1; %% este calculo funciona para 2 nod por elem
for inod = 1:Nnodporelem
    for idof = 1:Ndofpornod
        elemDof(e,idof+Ndofpornod*(inod-1)) = e*idof+(e-1)*Ndofporelem+(inod-1)*Ndofpornod;
    end
end
for e=2:N
    elemDof(e,1:Ndofpornod) = elemDof(e-1,Ndofpornod+1:end);
    elemDof(e,Ndofpornod+1:end) = elemDof(e,Ndofpornod)+1:elemDof(e,Ndofpornod)+Ndofpornod;
end

elemDof = elemDof +(nodstart-1)*Ndofpornod;
end

