%Placa Kirchoff
funcFormaKirch
N=shapefuns;

%% Problema - nodos/elementos

divisionesx = 6;
divisionesy = 6;
a=1;b=1; % Tama�o del problema
dx = a/(divisionesx-1);
dy = b/(divisionesy-1);

[nodos, elementos] = rectmesh(a,b,divisionesx,divisionesy);
%% DOFINITIONS
[Nelem, Nnodporelem]= size(elementos);  
[Nnod, Ndim] = size(nodos); % Numero de nodos, Numero de dimensiones del problema (es 2-D)

Ndofpornod = 3; %Para placa Kirchoff
dof = Nnod*Ndofpornod;
DOF = (1:dof)'; %vector columna

n2d = @(nodo) [nodo*3-2, nodo*3-1, nodo*3]; % Funci�n Node a DOF. Obtiene indices de dof de un nodo. Si hay mas/menos de 3 dof por nodo entonces cambia 
%% Bonus: Armo matriz elemDof para ensamblar Matriz rigidez rapido
elemDof = zeros(Nelem,Ndofpornod*Nnodporelem);
for e = 1:Nelem
   for n = 1:Nnodporelem
       elemDof(e,n2d(n)) = n2d(elementos(e,n));
   end
end

%% Condiciones de Borde y Cargas
isFixed = false(dof,1);
% simplementeApoyadoNodos = [1 2 3 4 5 10 15 20 25 24 23 22 21 16 11 6]';
for n =1:Nnod
    if nodos(n,1)==0 || nodos(n,2)==0 || nodos(n,2)==b || nodos(n,1)==a
        isFixed(n2d(n))=true;
    end
end
isFree = ~isFixed;

syms x y real

%% Propiedades del material
E = 1.000;
NU = 0.300;
t = 1.0;
F = E*t^3/(12*(1 - NU^2)); % Rigidez a la flexion.

C = F*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Obtencion Matriz Rigidez por integracion simbolica No-Isoparametrica
Kg = sparse(dof,dof);
Ke = double(int(int(B'*C*B,x,0,a/(divisionesx-1)),y,0,b/(divisionesy-1))); %Los elementos son 1x1. Integro asi
for e = 1:Nelem
    storeTo = false(dof,1);
    storeTo(elemDof(e,:)) = true;
    Kg(storeTo,storeTo) = Kg(storeTo,storeTo) + Ke;
end

%% Obtencion vector columna de cargas R
p0=1;
R = zeros(dof,1);
Nint = int(int(N.',x,0,dx),y,0,dy);
Q = double(p0*Nint);
for e = 1:Nelem 
    storeTo = false(dof,1);
    storeTo(elemDof(e,:))=true;
    R(storeTo)=R(storeTo)+Q;
end

%% Inversion (C�lculo de desplazamientos)
Dr = Kg(isFree,isFree)\R(isFree);


D=zeros(dof,1);
D(isFree) = Dr;

%% Grafico de Desplazamientos
Dz = zeros(divisionesx,divisionesy); % Matriz superficie
xv =[];
yv = [];
for n=1:Nnod
    xv = [xv nodos(n,1)];
    yv = [yv nodos(n,2)];
    Dz(n) = D(n*3-2);
end
Dz=reshape(Dz,[],1)';
scatter3(xv,yv,Dz)

