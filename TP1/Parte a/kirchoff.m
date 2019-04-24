%Placa Kirchoff

%% Problema - nodos/elementos
DIVITERk = ceil(1.2.^(1:14))+1

errveck = nan(length(DIVITERk),1);
Neleck = errveck;
timerk = Neleck;
iter = 1;
for divy = DIVITERk(2:end)
    
    iter = iter+1;
    tic
% divisionesx = 11; % Minimo 3 divisiones
% funcFormaMind8
if divy == DIVITERk(iter-1)
    continue
end
if mod(divy,2)==0
    divy=divy+1;

end
divisionesx = ceil(divy*1.4);
if mod(divisionesx,2)==0
    divisionesx=divisionesx+1;
end
% divisionesx = 11;
divisionesy = divy;
a=1.4;b=1; % Tamaño del problema
dx = a/(divisionesx-1);
dy = b/(divisionesy-1);

[nodos, elementos] = rectmesh(a,b,divisionesx,divisionesy);

funcFormaKirch
%% DOFINITIONS
[Nelem, Nnodporelem]= size(elementos);  
[Nnod, Ndim] = size(nodos); % Numero de nodos, Numero de dimensiones del problema (es 2-D)

Ndofpornod = 3; %Para placa Kirchoff
dof = Nnod*Ndofpornod;
DOF = (1:dof)'; %vector columna

n2d = @(nodo) [nodo*3-2, nodo*3-1, nodo*3]; % Función Node a DOF. Obtiene indices de dof de un nodo. Si hay mas/menos de 3 dof por nodo entonces cambia 
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
        isFixed(n2d(n))=[true false false];
    end
end
isFree = ~isFixed;

syms x y real

%% Propiedades del material
E = 210e9;
NU = 0.300;
t = a/100;
F = E*t^3/(12*(1 - NU^2)); % Rigidez a la flexion.

C = F*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Obtencion Matriz Rigidez por integracion simbolica No-Isoparametrica
Kg = sparse(dof,dof);
Ke = double(int(int(B'*C*B,x,[-dx dx]),y,[-dy dy]));



% jac = Nder*nodesEle;
% dNxy = jac\Nder;   % dNxy = inv(jac)*dN   
for e = 1:Nelem
    storeTo = elemDof(e,:);
    Ke = double(int(int(B'*C*B,x,[-dx dx]),y,[-dy dy]));
%     storeTo(elemDof(e,:)) = true;
    Kg(storeTo,storeTo) = Kg(storeTo,storeTo) + Ke;
end

%% Obtencion vector columna de cargas R
p0=-0.05e6; %Pa
R = zeros(dof,1);
Nint = int(int(N.',x,[-dx dx]),y,[-dy dy]);
Q = double(p0*Nint);
for e = 1:Nelem 
    storeTo = elemDof(e,:);
%     storeTo(elemDof(e,:))=true;
    R(storeTo)=R(storeTo)+Q;
end

%% Inversion (Cálculo de desplazamientos)
Dr = Kg(isFree,isFree)\R(isFree);
timerk(iter)=toc;

D=zeros(dof,1);
D(isFree) = Dr;

errveck(iter)=max(max(abs(D_err)));
% err =max(max(abs(D_err)));
% scatter(divy,err)
% hold on
% clear
Neleck(iter)=Nelem;

end
semilogx(Neleck,errveck)
title('Convergencia de solución')
ylabel('Error absoluto máximo [mm]')
xlabel('Numero de Elementos')

return
W = D(1:3:end);
fprintf("w_max = %f",max(abs(W)))
W_analytic = zeros(Nnod,1);
D_err = nan(Nnod,1);
N=9; %Iteraciones de sol analitica 9 es buen número, converge bastante bien
for n = 1:Nnod
    W_analytic(n) = w_analytic(nodos(n,1),nodos(n,2),a,b,N,p0,F);
    if ~(nodos(n,2)==0 || nodos(n,2)==b ||nodos(n,1)==0 ||nodos(n,1)==a)
        D_err(n) =  W_analytic(n) - W(n);
    end
end
figure
scatter3(nodos(:,1),nodos(:,2),D_err,'r')
errRel = max(abs(D_err))/max(abs(W_analytic));
title(sprintf("Error para placa espesor $t=%0.0f$mm Error relativo: %f\\%%",t*1000,errRel),'interpreter','latex')%Todo ese lio para imprimir un porcentaje
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
zlabel('$w$ [m]','interpreter','latex')
return
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

