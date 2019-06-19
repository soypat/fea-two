prog_time=tic;
oldpath=path;
path('auxfunc',oldpath);
path('meshed',path);
% HEXA 20 Shapefuns según ADINA
clear N dN dNaux NL ksi eta zeta x y z
clc

% aux = load('searelemCoarse.txt'); % ELEMENTOS
% elementos = aux(:,2:21);
% aux = load('searnodCoarse.txt');% NODOS
% nodos = aux(:,2:4);

aux = load('elemNX.txt'); % ELEMENTOS
elementos = aux(:,2:21); %H20
aux = load('nodNX.txt');% NODOS
nodos = aux(:,2:4);
nodenumbering=aux(:,1);
[elementos, nodos] = fix_it_ftlog(elementos,nodos,nodenumbering);


%% Funciones de Forma
funcforma

%% SEAR GEOMETRY
centro_agujero=[16-13.14 3/2 -2.425];
xnod_carga=16;


%% DOFINITIONS
Ndofpornod = 3; % x y z
[Nnod, Ndim]     =     size(nodos);
[Nelem, Nnodporelem] = size(elementos);
Ndofporelem  =  Ndofpornod*Nnodporelem;
dof = Ndofpornod*Nnod;
DOF = reshape(uint16(1:dof),Ndofpornod,[])'; %PARA DOF MAYOR A 65 MIL USAR DOUBLE
%% Material
E=200e3;
nu=.3;
%% Relacion constitutiva isótropa elastica
lambda = E*nu/(1+nu)/(1-2*nu);
G = E/(2+2*nu);
C = [lambda+2*G lambda lambda 0 0 0; %Relacion constitutiva 3D
     lambda lambda+2*G lambda 0 0 0;
     lambda lambda lambda+2*G 0 0 0;
        0     0     0     G    0  0;
        0     0     0     0    G  0;
        0     0     0     0    0  G];
 
    f=2+2*nu;
Cinv = [1 -nu -nu 0 0 0;
    -nu 1 -nu 0 0 0;
    -nu -nu 1 0 0 0;
    0 0 0 f 0 0;
    0 0 0 0 f 0;
    0 0 0 0 0 f]/E;
% Cinv == C^-1

%% Punto Gauss con método acelerado de calculo
[wpg, upg, npg] = gauss([3 3 3]); % Ver tabla 6.8-1 Cook
[~,~,dNauxs,~] = shapefunGP(upg,N,dN,dNaux,NL);

%% Matriz Rigidez Global

K = zeros(dof,dof);
fprintf('\nComienzo de calculos luego de %0.2f segundos\n\n',toc(prog_time))
k_time=tic;
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
%% Condicionar matriz a forza
% colSumsK=sum(K);
% isNullCol=colSumsK==0
tiempoK=toc(k_time);
digi=NaN;
cond_time=tic;
K=defuzz(K);
% digi=get_cond(K);
digi=1;

fprintf('\nMatriz de rigidez calculada en %0.2f segundos\nDigitos perdidos por condicionamiento: %0.1f\nTiempo en calcular condicionamiento: %0.1f\n',tiempoK,digi,toc(cond_time));
K=sparse(K);
qx = 10*33.5442; %33.5 MPa para alcanzar los 65 newtons sobre el area que es 1.9377mm^2
cargas

%% FIXITY
i=0;
fixed=false(dof,1);
k=0;
for inod = nodos' %Busco nodos en superficie
    x=inod(1);y=inod(2);z=inod(3);r=sqrt((x-centro_agujero(1))^2 + (z-centro_agujero(3))^2 );tita=atan2(z-centro_agujero(3),x-centro_agujero(1));
    i=i+1;
    if abs(r-2.75/2)<3e-2 %Condición z==0 (evitando error)
         fixed([i*3-2 i*3-1 i*3])=true;
         k=k+1;
    end
end

free=~fixed;
fprintf('\nComienzo de inversion\n\n')
% prepinv
% loose=abs(sum(K)).'==0;%Borra posiciones 'sueltas'
% free=and(free,~loose);
Rr=R(free);

% clear R
% Kr=K(free,free);


pause(.1);
inv_time = tic;

Dr = K(free,free)\Rr;

fprintf('\nDesplazamientos calculados en %0.0f segundos\n\n',toc(inv_time))
D = zeros(dof,1);
D(free)=Dr;

%% Recuperación de tensiones
S = zeros(Nelem,Nnodporelem,6);
strain = zeros(Nelem,Nnodporelem,6);
% unod declarado arriba!

%Funciones de formas evaluadas en los nodos
[Ns,~,dNauxs,~] = shapefunGP(unod,N,dN,dNaux,NL);
tension_time=tic;
for e = 1:Nelem
    index = elementos(e,:);
    elenod = nodos(index,:);
    for inod = 1:Nnodporelem
        % Derivadas de x,y, respecto de ksi, eta
        J    =  dNauxs{inod}*elenod;
        dNxyz = J\dNauxs{inod}; 
        
        B=zeros(size(C,2),size(Ke,1)); % B es 6x60
        B(1,1:Ndofpornod:end-2)=dNxyz(1,:); %  dNxy es 6x20 OJO
        B(2,2:Ndofpornod:end-1)=dNxyz(2,:); %  
        B(3,3:Ndofpornod:end)=dNxyz(3,:); %  dz %VER PAGINA 80 Cook. Ec (3.1-9)
        B(4,1:Ndofpornod:end-2)=dNxyz(2,:); %dy
        B(4,2:Ndofpornod:end-1)=dNxyz(1,:); %dx
        B(5,2:Ndofpornod:end-1)=dNxyz(3,:); %dz
        B(5,3:Ndofpornod:end)=dNxyz(2,:); %dy
        B(6,1:Ndofpornod:end-2)=dNxyz(3,:); %dz
        B(6,3:Ndofpornod:end)=dNxyz(1,:); %dx
        
        meindof = reshape(DOF(index,:)',1,[]);
        
        deformation = B*D(meindof);
        strain(e,inod,:) = deformation;
        S(e,inod,:) = C*deformation;% Stress
    end
end

fprintf('\nTensiones recuperadas en %0.2f segundos\n\n',toc(tension_time))
Svm=dametuvonvonFE(S);
figure(1)
bandplot3D(elementos,nodos,Svm,[],'k');

% save('results','elementos','nodos','S','D')
% clear
% res=open ('results.mat');D=res.D;elementos=res.elementos;nodos=res.nodos;S=res.S;
% figure(1)
% bandplot3D(elementos,nodos,S(:,:,1),[],'k');
figure(2)
escala=1;

posiciondeformada=nodos+escala*reshape(D,3,[]).';

bandplot3Dopaco(elementos,posiciondeformada,Svm,[],'k');
fprintf('\nPrograma finalizado luego de %0.2f segundos\n\n',toc(prog_time))
