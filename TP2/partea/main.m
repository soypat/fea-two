%% PROBLEM
LargoMotor = 5.719;
% Modifiables
AB = [1.8 1.8];
longE = [.7 .99 .8];%[1.2 2.95 ] %E1, E2, E3
rpmexc = 600;
omegaexc = rpmexc/60*2*pi;
Lelemax=5; %Longitud máxima de elementos
h = .09;
b = .04;

% No-modificables
rojoX = [417 667 1197 1727 2257 2787 3317 3847 4377 4907 5437]'/1000;

ancho = 1.93;
CG = [1.855 ancho/2 1];

%% Mesheador del basamento
meshMe
close all
% Draw_Barra(elementos,nodos,'k')
% title('Basamento')
%% Mesheo masa
nodos = [nodos;CG];
Nnod = size(nodos,1);
masa=Nnod;

Nvigas = size(elementos,1);
%% Vigamaterials
E = 200e9;
rho = 7900;
nu=0.3;

A=b*h;
Iz = b*h^3/12;
Iy = h*b^3/12;
% Jxy = h*b*(h^2 + b^2)/12
% Jp = h*b*h*b/64
Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
%% Dofinitions
Nnod = size(nodos,1);
n2d6=@(n) [n.*6-5 n.*6-4 n.*6-3 n.*6-2 n.*6-1 n.*6];
elemDof = [n2d6(elementos(:,1)) n2d6(elementos(:,2))];
dof = Nnod*6;
%% Mass Stiffness Matrices
K = zeros(dof,dof);
M = zeros(dof,dof);

for e = 1:Nvigas
    storeTo = elemDof(e,:);
    n1=nodos(elementos(e,1),:);n2=nodos(elementos(e,2),:);
    Le=norm(n2-n1);
    [ke, me] = vigastiffness(E,nu,rho,A,Iz,Iy,Jtors,Le);
    Ke = vigorotar(ke,n1,n2,[0 0 1]);
    Me = vigorotar(me,n1,n2,[0 0 1]);
    
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
    M(storeTo,storeTo) = M(storeTo,storeTo) + Me;
end

%% Reemplazo RigidLinks por Rigid Beams
isFixed = false(dof,1);
Nrigid = Nrojo*2;
krigid = 1e8;
rigidElementos =[];
[ke,~] = vigastiffness(E*100,nu,rho,A,Iz,Iy,Jtors,.5);
for e = 1:Nrojo %Unimos barras a la masa, para sustituir lo que seria un rigid link
    for i=1:2
    storeTo = [n2d6(rojos(e)) n2d6(masa)];
    n1=nodos(rojos(e,i),:);n2=nodos(masa,:);
    rigidElementos = [rigidElementos;rojos(e,i) masa];
    Ke = vigorotar(ke,n1,n2,[0 0 -1]);
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
%     isFixed(storeTo)= isFixed(storeTo) | [0 0 0 0 0 0 0 0 0 1 1 1]'; %Matamos giros en la masa puntual repetidas veces para asegurar su muerte
    end
end

%% ""Rigid Links""
% isFixed = false(dof,1);
% Nrigid = Nrojo*2;
% krigid = 1e8;
% rigidElementos =[];
% for e = 1:Nrojo %Unimos barras a la masa, para sustituir lo que seria un rigid link
%     for i=1:2
%     storeTo = [n2d6(rojos(e)) n2d6(masa)];
%     n1=nodos(rojos(e,i),:);n2=nodos(masa,:);
%     rigidElementos = [rigidElementos;rojos(e,i) masa];
%     [ke, ~] = vigastiffness(krigid,nu,0,1,0,0,0,1);
%     Ke = vigorotar(ke,n1,n2,[0 0 -1]);
%     K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
%     isFixed(storeTo)= isFixed(storeTo) | [0 0 0 0 0 0 0 0 0 1 1 1]'; %Matamos giros en la masa puntual repetidas veces para asegurar su muerte
%     end
% end
% Draw_Barra(rigidElementos,nodos(:,[1 3 2]),'k')
% title('Rigid Links plano xz')
% Draw_Barra(rigidElementos,nodos(:,[1 2 3]),'k')
% title('Rigid Links plano xy')
%% Acoplo masa puntual
masapuntual = 1e6*blkdiag(.055, .055, .055, 0.49, 0.28455,0.28455);
storeTo = n2d6(masa);
M(storeTo,storeTo)=M(storeTo,storeTo) + masapuntual;

% digitosPerdidos = get_cond(K);
% if digitosPerdidos >10
%     error('Condicionamiento malo!')
% end


%% Creo las condiciones de bordes en forma de bulones
Eb=200e9;
d =16e-3; %16mm
Ab= d^2/4*pi;
Lb =70e-3; %70mm
% Lbtransversal = d;
kb = Eb*Ab/Lb;
kbt = kb/10;
K = abulonar(azules(:,1),K,[kb kbt kbt]); %La restringo (como que aplico condiciones de borde)
K = abulonar(azules(:,2),K,[kb kbt kbt]); %Restringo las otras
%% Aplico condiciones de borde! En este caso es impedir que giren las barras conectadas a la masa puntual
isFree=~isFixed;
Kr = K(isFree,isFree);
Mr = M(isFree,isFree);
%% Busco Fuerzas para causar desplazamientos iniciales
dmasa = 10e-3;
% DEMOSTRACION:
% a = d^2x/dt^2 -> x = dmasa*sin(omegaexc*t) -> v = dmasa*omegaexc*cos(o*t)
% a = - dmasa*omegaexc^2 * sin(omegaexc*t) -> Fmax = dmasa*omegaexc^2 *masa
F0 = dmasa* omegaexc^2 * masapuntual(1,1);
Fmasa = [0 0 F0 0 0 0]; % se mueve 10mm en z
storeTo = n2d6(masa);
R = zeros(dof,1);
R(storeTo) = Fmasa;
Rr=R(isFree);
%% Cuasiestatico
Fcuasi = [0 0 -masapuntual(1,1)*9.81 0 0 0];
Rcuasi = zeros(dof,1);
Rcuasi(n2d6(masa)) = Fcuasi;
Dr = Kr\Rcuasi(isFree);
D = zeros(dof,1);
D(isFree)=Dr;
maxsig=0;

for e = 1:Nvigas
    
    n1 = nodos(elementos(e,1),:);
    n2 = nodos(elementos(e,2),:);
    T = Tv(n2-n1,[0 0 1]);
    [ke, me] = vigastiffness(E,nu,rho,A,Iz,Iy,Jtors,Le);
    Ke =T*ke*T';
    Dlocal=D(elemDof(e,:));
    flocal=ke*T'*Dlocal;
    Nx =  flocal(1);
    Vy =  flocal(2);
    Vz =  flocal(3);
    Mx1 = flocal(4);
    My1 = flocal(5);
    Mz1 = flocal(6);
    Mx2 = flocal(10);
    My2 = flocal(11);
    Mz2 = flocal(12);
    Mz=max([Mz1 Mz2]);
    My = max([My1 My2]);
    A=b*h;
    sig = Nx/A+abs(Mz*h/(2*Iz));
    sig2 = Nx/A+abs(My*b/(2*Iy));
    if maxsig<max([sig,sig2]) && max([sig,sig2])<1e9
        maxsig=max([sig,sig2]);
        badel=e;
    end
end


% resonance





