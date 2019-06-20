clear N dN dNaux NL ksi eta zeta x y z
clc
%%PROBLEM
Q = 2e3; %2kW/m^3 generados
L=0.8;

%% SOLVE:
funcforma

div = 3;
[nodos, elementos] = mesh3D([0 L;0 L;0 L],[div div div]);
scatter3(nodos(:,1),nodos(:,2),nodos(:,3))
medio = ceil(div^3/2);

Ndofpornod = 1; % T
[Nnod, Ndim]     =     size(nodos);
[Nelem, Nnodporelem] = size(elementos);
Ndofporelem  =  Ndofpornod*Nnodporelem;
dof = Ndofpornod*Nnod;
DOF = reshape(1:dof,Ndofpornod,[])';

%% Material
cp=528; %J/kg.K
rho =4500;
k=17;%W/m/K
C = [k 0 0 
     0 k 0
     0 0 k]; %considero material isotropico

%% Punto Gauss con método acelerado de calculo
[wpg, upg, npg] = gauss([2 2 2]); % Ver tabla 6.8-1 Cook
[Ns,~,dNauxs,~] = shapefunGP(upg,N,dN,dNaux,NL);

%% Rigidity Assemble and charges
K = zeros(dof,dof);
R=zeros(dof,1);
for e = 1:Nelem
    Ke = zeros(Ndofporelem);
    r=zeros(Ndofporelem,1);
    index = elementos(e,:);
    elenod = nodos(index,:);
    for ipg = 1:npg
        J    = dNauxs{ipg}*elenod;
        dNxyz = J\dNauxs{ipg};
        B = [dNxyz(1,:);dNxyz(2,:);dNxyz(3,:)];
        
        Ke = Ke + B'*C*B*wpg(ipg)*det(J);
        r = r+Ns{ipg}'*Q*wpg(ipg)*det(J);
    end
    meindof = reshape(DOF(index,:)',1,[]);
    R(meindof) = R(meindof)+r; %cargas térmicas por generacion
    
    K(meindof,meindof)=K(meindof,meindof)+Ke;
end

T = zeros(dof,1);
T(medio) = 300; %K condicion de borde

cc = false(dof,1);
cc(medio)=true;
xx = ~cc;

Tx = K(xx,xx)\(R(cc)-K(xx,cc)*T(cc));
Qx = K(cc,xx)*T(xx)+K(cc,cc)*T(cc)
R(cc)=Qx
