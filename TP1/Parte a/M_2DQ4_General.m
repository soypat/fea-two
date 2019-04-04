% Modelo Q4 Tracci�n simple. Integraci�n por Gauss.
% clear
close all
format short g

% Discretizacion
nodes = [ 0.0  0.0
          1.0  0.0
          2.0  0.0
          0.0  1.5
          1.25 0.75
          2.0  1.0
          0.0  2.0
          1.0  2.0
          2.0  2.0 ];

% nodes = [ 0.0  0.0
%           1.0  0.0
%           2.0  0.0
%           0.0  1.0
%           1.0  1.0
%           2.0  1.0
%           0.0  2.0
%           1.0  2.0
%           2.0  2.0 ];      
      
elements = [1  2  5  4
            2  3  6  5
            4  5  8  7
            5  6  9  8];      %Matriz de conectividades

nDofNod = 2;                    % N�mero de grados de libertad por nodo
nNodEle = 4;                    % N�mero de nodos por elemento
nel = size(elements,1);         % N�mero de elementos
nNod = size(nodes,1);           % N�mero de nodos
nDofTot = nDofNod*nNod;         % N�mero de grados de libertad

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc(1,1:2) = true;
bc(2,2) = true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R([7 9],2) = 0.5;
R(8,2) = 1.0;
R(3,2) = -0.5;

% Propiedades del material
E = 1;
NU = 0.3;

meshplot(elements,nodes,'b')

%% Matriz Constitutiva (plane stress)

C = E/(1 - NU^2)*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];
               

%% Gauss           
a   = 1/sqrt(3);
% Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
% N�mero de puntos de Gauss
npg = size(upg,1);
wpg = ones(npg,1);

%% Matriz de rigidez
K = zeros(nDofTot);
nodeDofs = reshape(1:nDofTot,nDofNod,nNod)';
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodesEle = nodes(elements(iele,:),:);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:); 

        Ke = Ke + B'*C*B*wpg(ipg)*det(jac);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Reduccion Matriz
isFixed = reshape(bc',[],1);
isFree = ~isFixed;

Rr = reshape(R',[],1);

% Solver
Dr = K(isFree,isFree)\Rr(isFree);

% Reconstrucci�n
D = zeros(nDofTot,1);
D(isFree) = D(isFree) + Dr;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones = nan(nDofTot,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,nDofNod,[]))';

%% Recuperaci�n de tensiones en los nodos
stress = zeros(nel,nNodEle,3);
uNod = [ -1 -1
          1 -1
          1  1
         -1  1 ];
for iele = 1:nel
    nodesEle = nodes(elements(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
                  -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:2:7) = dNxy(1,:);
        B(2,2:2:8) = dNxy(2,:);
        B(3,1:2:7) = dNxy(2,:);
        B(3,2:2:8) = dNxy(1,:);
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*B*D(eleDofs);
    end
end

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficaci�n
bandplot(elements,nodePosition,stress(:,:,2),[],'k');
meshplot(elements,nodes,'b')
