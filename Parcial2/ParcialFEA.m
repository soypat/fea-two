clc
clearvars

%% malla

nodes = [ 0 0
          1 0
          2 0
          0 1
          1 1
          2 1
          0 2
          1 2
          2 2];
      
elements  = [ 1 2 5 4
              2 3 6 5
              4 5 8 7
              5 6 9 8];
          
eleType = 'Q4';

meshplot(elements,nodes,nodes*0,'b',1);

%% parametros del modelo

nDofNod = 1;                    % N?mero de grados de libertad por nodo
nel = size(elements,1);         % N?mero de elementos
nNod = size(nodes,1);           % N?mero de nodos
nEle = size(elements,1);
nNodEle = size(elements,2);     % N?mero de nodos por elemento
nDofTot = nDofNod*nNod;         % N?mero de grados de libertad


%% Conductividad
Kappa = 0.8*eye(2);

%% espesor
t = 1; 
%% Capacidad
% C = eye(2);

%% Pungos de Gauss

% Gauss

a   = 1/sqrt(3);
        % Ubicaciones puntos de Gauss
upg = [ -a  -a
         a  -a
         a   a
        -a   a ];    
    
npg = size(upg,1);
wpg = ones(npg,1);

dof = reshape(1:nNod*nDofNod, nDofNod, nNod)';

K = zeros(nDofTot);
C = zeros(nDofTot);

for iele = 1:nEle
    nodesEle = nodes(elements(iele,:),:);
    eleDofs = reshape(dof(elements(iele,:),:)', [], 1);
    Ke = zeros(nDofNod*nNodEle);
    for ipg = 1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);  
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta], eleType);  
        N = shapefuns([ksi eta], eleType);

        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;                  
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        B = zeros(2,4);
        B(1,:) = dNxy(1,:);
        B(2,:) = dNxy(2,:);
%         Ce = Ce + B'*C_ele*B*wpg(ipg)*t*det(jac);
        Ke = Ke + B'*Kappa*B*wpg(ipg)*t*det(jac);
    end
    K(eleDofs, eleDofs) = K(eleDofs, eleDofs) + Ke;
%     C(eleDofs, eleDofs) = C(eleDofs, eleDofs) + Ce;
end

%Condiciones de Borde
Q=[0 0 0 0 0 0 0.5 1 0.5]';

bc=[1 1 1 0 0 0 0 0 0];

isfree=~bc;

T=zeros(nDofTot,1)

%Resolución
T(isfree)=K(isfree,isfree)\Q(isfree);

%Reorganizacion resultados
for iele=1:nel
    Tel(iele,:)=T(elements(iele,:));
end

%Ploteo
bandplot(elements,nodes,Tel,[min(T) max(T)],'k',16)