funcFormaMind8

%% Problema - nodos/elementos

divisionesx = 25; % Minimo 3 divisiones
divisionesy = 17; % Minimo 3 divisiones
a=1.400; b=1.000; % Tamaño del problema
dx = a/(divisionesx-1);
dy = b/(divisionesy-1);

[nodos, elementos] = rectmeshQ8(a,b,divisionesx,divisionesy);

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
E = 210e9; %GPa Acero
NU = 0.3;
t = a/100; % a mm

F = E*t^3/(12*(1 - NU^2)); %Rigidez ante la flexion
G = E/(2+2*NU); % Rigidez a la torsion
Cb = -  [F NU*F 0;
    NU*F F 0 
    0 0 (1-NU)*F/2];
Cs = -5/6*[G*t 0;0 G*t];
C = blkdiag(Cb,Cs);

%% Gauss  2x2      
k   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg2 = [ -k  -k
         k  -k
         k   k
        -k   k ];
npg2 = size(upg2,1);
wpg2 = ones(npg2,1);
%% Gauss 3x3 - Full integration
k=sqrt(0.6);
upg3 = [-k -k;-k 0; -k k;0 -k;0 0;0 k;k -k;k 0;k k];
k = [5/9;8/9;5/9];
wpg3 = reshape(k*k',1,[]);
npg3 = length(wpg3);
%% Allocate N, B and dN. (va mucho mas rapido con esto)
% GAUSS 3x3
Ns3 = cell(npg3,1);
dNs3 = cell(npg3,1);
for ipg = 1:npg3
        ksi = upg3(ipg,1); eta = upg3(ipg,2);
        dNs3{ipg} = double(subs(dN));
        Ns3{ipg} = double(subs(N));
end
Ns2 = cell(npg2,1);
dNs2 = cell(npg2,1);
for ipg = 1:npg2
        ksi = upg2(ipg,1); eta = upg2(ipg,2);
        dNs2{ipg} = double(subs(dN));
        Ns2{ipg} = double(subs(N));
end
%% Obtencion Matriz Rigidez por integracion Isoparametrica
Kg = sparse(dof,dof);
for e = 1:Nelem
    Kb = zeros(Ndofporelem);
    Ks = zeros(Ndofporelem);
    storeTo = false(dof,1);
    storeTo(elemDof(e,:)) = true;
    nodesEle = nodos(elementos(e,:),:);
    Bb = zeros(3,Ndofporelem);
    Bs = zeros(2,Ndofporelem);
    %% K Bending
    for ipg = 1:npg3
        ksi = upg3(ipg,1); eta = upg3(ipg,2);

        Nder=dNs3{ipg};
        
        jac = Nder*nodesEle;
        dNxy = jac\Nder;   % dNxy = inv(jac)*dN     
        for i = 1:Nnodporelem % Armo matriz B de bending
            Bb(:,(i*3-2):(i*3)) = [0 dNxy(1,i) 0
                0 0 dNxy(2,i)
                0 dNxy(2,i) dNxy(1,i)];
        end
        Kb = Kb + Bb'*Cb*Bb*wpg3(ipg)*det(jac);
    end
    %% K shear
    for ipg = 1:npg2
    ksi = upg2(ipg,1); eta = upg2(ipg,2);

        Nder=dNs2{ipg};
        
        jac = Nder*nodesEle;
        dNxy = jac\Nder;   % dNxy = inv(jac)*dN
        
        for i = 1:Nnodporelem 
            Bs(:,(i*3-2):(i*3)) = [-dNxy(1,i) Ns2{ipg}(i) 0
                -dNxy(2,i) 0 Ns2{ipg}(i)]; % Ver cook 15.3-3
        end
        
        Ks = Ks + Bs'*Cs*Bs*wpg2(ipg)*det(jac);
    end
    %% Acople general
    Kg(storeTo,storeTo) = Kg(storeTo,storeTo) + Kb + Ks;
end

%% Condiciones de Borde (empotrado)
isFixed = false(dof,1);
for n = 1:Nnod
   x = nodos(n,1);y = nodos(n,2);
   if x==0 || y ==0 || y == b || x ==a
       isFixed(n2d(n))=[true false false];
   end
end
isFree = ~isFixed;

%% Cargas
p0 = -0.071e6 ; %Pa

R = zeros(dof,1);

for e = 1:Nelem
    storeTo = false(dof,1);
    storeTo(elemDof(e,1:3:end))=true;
    nodesEle = nodos(elementos(e,:),:);
    for ipg = 1:npg
        jac = dNs{ipg}*nodesEle;
        Q = Ns{ipg}'*p0*wpg(ipg)*det(jac);
        R(storeTo)=R(storeTo)+Q;
    end
end
%% Solucion
Dr = Kg(isFree,isFree)\R(isFree);


D=zeros(dof,1);
D(isFree) = Dr;

%% Graficar
placaplotQ8(nodos,D)
% Dz = zeros(divisionesx,divisionesy); % Matriz superficie
% xv =[];
% yv = [];
% for n=1:Nnod
%     xv = [xv nodos(n,1)];
%     yv = [yv nodos(n,2)];
%     Dz(n) = D(n*3-2);
% end
% Dz=reshape(Dz,[],1)';
% scatter3(xv,yv,Dz)
% hold on
% % figure 
% xv = 0:dx:a;
% yv = 0:dy:b;
% zv = reshape(Dz,length(xv),length(yv));
% [X, Y] = meshgrid(xv,yv);
% surf(X,Y,zv,'FaceAlpha',0.5)