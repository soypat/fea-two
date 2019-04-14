funcFormaMind

%% Problema - nodos/elementos

divisionesx = 11;
divisionesy = 7;
a=1.4; b=1; % Tamaño del problema
dx = a/(divisionesx-1);
dy = b/(divisionesy-1);

[nodos, elementos] = rectmesh(a,b,divisionesx,divisionesy);

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

%% Propiedades Material
E = 210e9; %GPa Aluminio
NU = 0.3;
t = a/100; % a mm

F = E*t^3/(12*(1 - NU^2)); %Rigidez ante la flexion
G = E/(2+2*NU); % Rigidez a la torsion
Cb = [F NU*F 0;
    NU*F F 0 
    0 0 (1-NU)*F/2];
Cs = 5/6*[G*t 0;0 G*t];
C = blkdiag(Cb,Cs);

%% Gauss 2x2       
k   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg2 = [ -k  -k
         k  -k
         k   k
        -k   k ];
npg2 = size(upg2,1);
wpg2 = ones(npg2,1);
%% Gauss 1x1
upg1 = [0 0];
wpg1 = 4;
npg1 = 1;

%% Allocate N, B and dN. (va mucho mas rapido con esto)
% Bs = cell(npg2,1);
Ns2 = cell(npg2,1);
dNs2 = cell(npg2,1);
for ipg = 1:npg2
        ksi = upg2(ipg,1); eta = upg2(ipg,2);
        dNs2{ipg} = double(subs(dN));
        Ns2{ipg} = double(subs(N));
end
% Allocate Gauss 1x1
ksi = upg1(1); eta = upg1(2);
dNs1 = double(subs(dN));
Ns1 = double(subs(N));

%% Obtencion Matriz Rigidez por integracion simbolica No-Isoparametrica
Kg = sparse(dof,dof);
for e = 1:Nelem
    Kb = zeros(Ndofporelem);
    dofIndex = elemDof(e,:);
    storeTo = elemDof(e,:);
    nodesEle = nodos(elementos(e,:),:);
    Bb = zeros(3,Ndofporelem);
    Bs = zeros(2,Ndofporelem);
    for ipg = 1:npg2
        Nder=dNs2{ipg};
        jac = Nder*nodesEle;
        dNxy = jac\Nder;   % dNxy = inv(jac)*dN
        for i = 1:Nnodporelem 
            Bb(:,(i*3-2):(i*3)) = [0 dNxy(1,i) 0;0 0 dNxy(2,i);0 dNxy(2,i) dNxy(1,i);];
        end
        Kb = Kb + Bb'*Cb*Bb*wpg2(ipg)*det(jac);
    end
    jac = dNs1*nodesEle;
    dNxy = jac\dNs1;   % dNxy = inv(jac)*dN
    for i = 1:Nnodporelem 
        Bs(:,(i*3-2):(i*3)) = [-dNxy(1,i) Ns1(i) 0;-dNxy(2,i) 0 Ns1(i)];
    end
    Ks = Bs'*Cs*Bs*wpg1*det(jac);

    Kg(storeTo,storeTo) = Kg(storeTo,storeTo) + Kb+ Ks;
end

%% Condiciones de Borde (empotrado)
isFixed = false(dof,1);
for n = 1:Nnod
   x = nodos(n,1);y = nodos(n,2);
   if  y ==0 || y == b || x==0 || x ==a 
       isFixed(n2d(n))=[true true true];
   end
end
isFree = ~isFixed;

%% Cargas
p0 = -0.05e6; %MPa

R = zeros(dof,1);

for e = 1:Nelem
    storeTo = elemDof(e,1:3:end);
%     storeTo(elemDof(e,1:3:end))=true;
    nodesEle = nodos(elementos(e,:),:);
    for ipg = 1:npg2
        jac = dNs2{ipg}*nodesEle;
        Q = Ns2{ipg}'*p0*wpg2(ipg)*det(jac);
        R(storeTo)=R(storeTo)+Q;
    end
end

%% Solucion
Dr = Kg(isFree,isFree)\R(isFree);


D=zeros(dof,1);
D(isFree) = Dr;

W= D(1:3:end);

fprintf("w_maxQ4 = %f",max(abs(W)))
%% Graficar
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
hold on
% figure 
% xv = 0:dx:a;
% yv = 0:dy:b;
% zv = reshape(Dz,length(xv),length(yv));
% [X, Y] = meshgrid(xv,yv);
% surf(X,Y,zv,'FaceAlpha',0.5)
%% Flo compare
% FloHead=open('flo.mat');
% Df=FloHead.Df;
% Dzf=FloHead.Dzf;
% Kf=FloHead.Kf;
% Rf = FloHead.Rf;
% libre=FloHead.libre;
% Sbf = FloHead.Sbf;
% 
% 
% % scatter3(xv,yv,Dzf)
% % legend('Patty','Flo')