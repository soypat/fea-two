function [nodos,elementos] = mesh3D(X,div)
%makes 3D mesh for H8 element - patty got u covered
% X = [x1 x2; y1 y2; z1 z2];
if sum(div<=1)>0
    error('ingrese division mayor a 1')
end
Nnod = prod(div);
nodos = zeros(Nnod,3);
gv = cell(3,1);
L = zeros(3,1);
dElem = zeros(3,1);
for i=1:3
    L(i) = X(i,2)-X(i,1);
    gv{i} = 0:L(i)/(div(i)-1):L(i);
    dElem(i) = div(i)-1;
end
n=0;
for nz = 1:div(3)
   for ny = 1:div(2)
       for nx = 1:div(1)
           n=n+1;
          xgv=gv{1}; ygv=gv{2}; zgv = gv{3};
          nodos(n,:)=[xgv(nx) ygv(ny) zgv(nz)];
       end
   end
end
Nnodpelem = 8;

Nelem = prod(dElem);
elementos = zeros(Nelem,Nnodpelem);
xyTot = div(1)*div(2);
e=1;
for ez = 1:dElem(3)
for ey = 1:dElem(2)
for ex = 1:dElem(1)
        firstnode = ex+div(1)*(ey-1)+(ez-1)*div(1)*div(2);
        znod = (xyTot+firstnode);
        index = [firstnode firstnode+1 firstnode+div(1)+1 firstnode+div(1) znod znod+1 znod+div(1)+1 znod+div(1)];
        elementos(e,:)=index;
        e=e+1;
end
end
end
scatter3(nodos(:,1),nodos(:,2),nodos(:,3))


end

