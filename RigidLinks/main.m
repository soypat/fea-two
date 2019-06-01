%% Rigid Link this One


nodos = [0 0 0
         1 0 0
         2 0 0
         3 0 0
         4 0 0
         1.3 1 1];
masa = 5;
empot= 3;
roller = 5;
rl = [2 4];

elementos = [1 2;2 3;3 4;4 5];
links = [2 6;4 6];

%% Dofinitions

Ndofpornod=6;
[Nelem, Nnodporelem] = size(elementos);
[Nnod, Ndim] = size(nodos);


elemDof=node2dof(elementos,Ndofpornod);


for e=1:Nelem
    
end

