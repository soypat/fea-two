% Make square FEM mesh
function [nodos, elementos] = rectmesh(a,b,divx,divy)
% divx = 4;
% divy = 3;
% a = 2;
% b = 3;

Nnod = divx*divy;

xgv = 0:a/(divx-1):a;
ygv = 0:b/(divy-1):b;
nodos = zeros(Nnod,2);
n=0;
for ny = 1:divy
   for nx = 1:divx
       n=n+1;
      nodos(n,:)=[xgv(nx) ygv(ny)];
       
   end
end
Nnodporelem = 4;
xelem = (divx-1);
yelem = (divy-1);
Nelem = xelem*yelem;

elementos = zeros(Nelem,Nnodporelem);
for ey = 1:yelem
    for ex = 1:xelem
        elemplace = ex+xelem*(ey-1);
        firstnode = ex+divx*(ey-1);
        index = [firstnode firstnode+1 firstnode+divx+1 firstnode+divx];
        elementos(elemplace,:) = index;
    end
end
