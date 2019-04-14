% Make square FEM mesh
function [nodos, elementos] = rectmeshQ8(a,b,divx,divy)
% RECTMESHQ8  rectmeshQ8(a,b,divx,divy)  div>=2
Nnodporelem = 8;

% Nnod = (2*divx-1)*(divy*2-1)-(divx-1)-(divy-1);

xgv = 0:a/(2*divx-2):a;
ygv = 0:b/(2*divy-2):b;
nodos = [];
n=0;
for ny = 1:length(ygv)
    for nx = 1:length(xgv)
        if mod(ny,2)==0 %Entonces cae sobre linea de centro Q8
            if mod(nx,2)==0 % Entonces estoy sobre el centro de un elemento. No hay nodo
                continue
            end
        end
        n=n+1;
        nodos=[nodos;xgv(nx) ygv(ny)];
    end
end
Nnod = size(nodos,1);

xelem = (divx-1);
yelem = (divy-1);
Nelem = xelem*yelem;

elementos = zeros(Nelem,Nnodporelem);
% uNod = [-1 -1
%         1 -1
%         1 1
%         -1 1
%         0 -1
%         1 0
%         0 1
%         -1 0];
myIndex = [1 3 5 7 2 4 6 8];
% ConvertToADINA = [1 3 5 7 2 4 6 8];
for ey = 1:yelem
    for ex = 1:xelem
        elemplace = ex+xelem*(ey-1);
        firstnode = 2*ex-1+(3*divx-1)*(ey-1);
        index = [firstnode firstnode+1 firstnode+2 firstnode+2*divx-ex+1 firstnode+3*divx+1 firstnode+3*divx firstnode+3*divx-1 firstnode+2*divx-ex];
        index = index(myIndex);
        elementos(elemplace,:) = index;
    end
% uNod == nodos(elementos,:)
% un=1
end