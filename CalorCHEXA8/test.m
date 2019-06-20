clear N dN dNaux NL ksi eta zeta x y z
clc

funcforma

nodos = unod*.4;
elementos = 1:8;

Ndofpornod = 1; % T
[Nnod, Ndim]     =     size(nodos);
[Nelem, Nnodporelem] = size(elementos);
Ndofporelem  =  Ndofpornod*Nnodporelem;
dof = Ndofpornod*Nnod;
DOF = reshape(1:dof,Ndofpornod,[])';

%% Material
E=200e3;
nu=.3;
k=30;%W/m/K
C = [k -k
     -k k];

%% Punto Gauss con método acelerado de calculo
[wpg, upg, npg] = gauss([2 2 2]); % Ver tabla 6.8-1 Cook
[~,~,dNauxs,~] = shapefunGP(upg,N,dN,dNaux,NL);

%% Rigidity Assemble
K = zeros(dof,dof);

for e = 1:Nelem
    Ke = zeros(Ndofporelem);
    index = elementos(e,:);
    elenod = nodos(index,:);
    for ipg = 1:npg
        J    = dNauxs{ipg}*elenod;
        dNxyz = J\dNauxs{ipg};
        
        B=zeros(size(C,2),size(Ke,1)); % B es 6x60
        B(1,1:Ndofpornod:end-2) = dNxyz(1,:); %  dNxy es 6x20 OJO
        B(2,2:Ndofpornod:end-1) = dNxyz(2,:); %  
        B(3,3:Ndofpornod:end)   = dNxyz(3,:); %  dz %VER PAGINA 80 Cook. Ec (3.1-9)
        B(4,1:Ndofpornod:end-2) = dNxyz(2,:); %dy
        B(4,2:Ndofpornod:end-1) = dNxyz(1,:); %dx
        B(5,2:Ndofpornod:end-1) = dNxyz(3,:); %dz
        B(5,3:Ndofpornod:end)   = dNxyz(2,:); %dy
        B(6,1:Ndofpornod:end-2) = dNxyz(3,:); %dz
        B(6,3:Ndofpornod:end)   = dNxyz(1,:); %dx

        Ke = Ke + B'*C*B*wpg(ipg)*det(J);
    end
    meindof = reshape(DOF(index,:)',1,[]);
    K(meindof,meindof)=K(meindof,meindof)+Ke;
end




%% Relacion constitutiva isótropa elastica
% lambda = E*nu/(1+nu)/(1-2*nu);
% G = E/(2+2*nu);
% C = [lambda+2*G lambda lambda 0 0 0; %Relacion constitutiva 3D
%      lambda lambda+2*G lambda 0 0 0;
%      lambda lambda lambda+2*G 0 0 0;
%         0     0     0     G    0  0;
%         0     0     0     0    G  0;
%         0     0     0     0    0  G];