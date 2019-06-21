clear N dN dNaux NL ksi eta zeta x y z
clc
%%PROBLEM
Q = 2e3; %2kW/m^3 generados
L=0.8;

%% SOLVE:
funcforma

div = 5;
[nodos, ~,elementos] = mesh3D([0 L;0 L;0 L],[div div div]);
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
wallnod=false(Nnod,1);
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
for n=1:Nnod
    x=nodos(n,1);y=nodos(n,2);z=nodos(n,3);
    if x==0||x==.8 || z==.8 || z==0 || y==0 || y==.8
        wallnod(n)=true;
    end
end


T = zeros(dof,1);

%% Aplicamos condiciones de borde conocidas 
%  para obtener un pantallazo del problema
cc = wallnod;
T(cc) = 2.7;
cc(medio)=true;
T(medio) = 300; %K condicion de borde
xx = ~cc;

Tx = K(xx,xx)\(R(xx)-K(xx,cc)*T(cc));
Qx = K(cc,xx)*Tx+K(cc,cc)*T(cc);
R(cc)=Qx;
T(xx)=Tx;

%% Iteracionar al radiacion
Rrad=zeros(dof,1);
Tespacio = 2.7;
Trad = Tespacio.^4;
boltz = 1.417e-8;
supnod = [1 2 3 4
          1 5 8 4
          4 8 7 3
          3 7 6 2
          2 6 5 1
          5 6 7 8];
k=1/sqrt(3);
wpgs=ones(4,1);
upg2d=[-k -k;-k k;k k;k -k];
upg1 = [-k -k k k]';
upg2 = [-k k k -k]';
upgS = ones(4,1);
npgs=4;
olds =0;
for e = 1:Nelem
    index = elementos(e,:);
    elenod = nodos(index,:);
    s=0;
    for snod = supnod'
        s=s+1;
        xnod=elenod(snod,1); ynod=elenod(snod,2); znod=elenod(snod,3);
        if sum(  xnod==0 | xnod==0.8 | ynod==0 | znod ==.8 | znod==0  )>2
            if olds~=s
            switch s
                case 1
                    surf_upg = [upg1 upg2 -upgS];
                case 2
                    surf_upg = [upg1 -upgS upg2];
                case 3
                    surf_upg = [upgS upg1  upg2];
                case 4
                    surf_upg = [upg1 upgS upg2];
                case 5
                    surf_upg = [-upgS upg1  upg2];
                case 6
                    surf_upg = [upg1 upg2 upgS];
            end
            [Ns, ~ ,dNauxs, ~ ] = shapefunGP(surf_upg,N,dN,dNaux,NL);
            end
            k=k+1;
            r = zeros(4,1);
            supindex = index(supnod(s,:));
            for ipg = 1:npgs %Acá comienza la integracion sobre superficie
                J = dNauxs{ipg}*elenod;
                
                Tupg = T(index)'*Ns{ipg}'; %interpolacion
                r = r-Ns{ipg}(supnod(s,:))'*boltz*(Tupg^4 - Trad)*det(J)*wpg(ipg);
                
            end
            Rrad(supindex)=Rrad(supindex)+r;
        end
    end
end



%% Graph problema

% scatter3(nodos(:,1),nodos(:,2),nodos(:,3))
% hold on
% scatter3(nodos(medio,1),nodos(medio,2),nodos(medio,3),'r*')
% hold on
% scatter3(nodos(wallnod,1),nodos(wallnod,2),nodos(wallnod,3),'m.')