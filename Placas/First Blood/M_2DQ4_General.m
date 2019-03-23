% Modelo Cascara Kirchoff Integraci�n por Gauss.
% clear
close all
format short 
clear dN N dNaux
funcForma

% Discretizacion
% nodes = [ 0.0  0.0
%           1.0  0.0
%           2.0  0.0
%           0.0  1.5
%           1.25 0.75
%           2.0  1.0
%           0.0  2.0
%           1.0  2.0
%           2.0  2.0 ];


nodes = [ 0.0  0.0
          1.0  0.0
          2.0  0.0
          0.0  1.0
          1.0  1.0
          2.0  1.0
          0.0  2.0
          1.0  2.0
          2.0  2.0 ];      
      
elements = [1  2  5  4
            2  3  6  5
            4  5  8  7
            5  6  9  8];      %Matriz de conectividades

nDofNod = 3;                    % N�mero de grados de libertad por nodo
nNodEle = 4;                % N�mero de nodos por elemento

nel = size(elements,1);         % N�mero de elementos
nNod = size(nodes,1);           % N�mero de nodos

nDofTot = nDofNod*nNod;         % N�mero de grados de libertad
DOF = 1:nDofTot;

bc = false(nNod,nDofNod);       % Matriz de condiciones de borde
bc([1 2 3 4 6 7 8 9],:)=true;

R = zeros(nNod,nDofNod);        % Vector de cargas
R(5,3) = .7081;


% Propiedades del material
E = 1;
NU = 0.3;
t = 1;
meshplot(elements,nodes,'b')

%% Matriz Constitutiva (plane stress)


C = E*t^3/12/(1 - NU^2)*[ 1.0     NU         0.0
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
        x = upg(ipg,1);
        y = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        Bs = double(subs(B));

        Ke = Ke + Bs'*C*Bs*wpg(ipg);
    end
    eleDofs = nodeDofs(elements(iele,:),:);
    eleDofs = reshape(eleDofs',[],1);
    K(eleDofs,eleDofs) = K(eleDofs,eleDofs) + Ke;  
end

%% Reduccion Matriz
% 
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
        x = uNod(inode,1);
        y = uNod(inode,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
%         dN = 1/4*[-(1-eta)   1-eta    1+eta  -(1+eta)
%                   -(1-ksi) -(1+ksi)   1+ksi    1-ksi ];  
        % Derivadas de x,y, respecto de ksi, eta
%         jac = dN*nodesEle;                      
        % Derivadas de las funciones de forma respecto de x,y.
%         dNxy = jac\dN;          % dNxy = inv(jac)*dN
        Bs = double(subs(B));
        
        eleDofs = nodeDofs(elements(iele,:),:);
        eleDofs = reshape(eleDofs',[],1);
        stress(iele,inode,:) = C*Bs*D(eleDofs);
    end
end

%% Configuracion deformada
D = (reshape(D,nDofNod,[]))';
nodePosition = nodes + D(:,1:2);

%Graficaci�n
bandplot(elements,nodePosition,stress(:,:,2),[],'k');
meshplot(elements,nodes,'b')
