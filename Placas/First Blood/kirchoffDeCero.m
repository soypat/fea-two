%Placa Kirchoff
funcForma
N=shapefuns;

%% Problema - nodos/elementos
nodes = [ 0.0  0.0
          1.0  0.0
          2.0  0.0
          0.0  1.0
          1.0  1.0
          2.0  1.0
          0.0  2.0
          1.0  2.0
          2.0  2.0 ]-1;      
      
elements = [1  2  5  4
            2  3  6  5
            4  5  8  7
            5  6  9  8];     
%% DOFINITIONS
[Nelem,Nnodporelem ]= size(elem);
[Nnod, Ndim] = size(nodes);

Ndofpornod = 3; %Para placa Kirchoff
dof = Nnod*Ndofpornod;
DOF = (1:dof)'; %vector columna

%% Condiciones de Borde y Cargas
isFree = ;

%% Propiedades del material
E = 1;
NU = 0.3;
t = 1;
C = E*t^3/12/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
               
 
