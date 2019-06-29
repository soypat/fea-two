%% PROBLEM
LargoMotor = 5.719;
% Modifiables
AB = [1.8 1.8];
longE = [.7 1 .8];%[1.2 2.95 ] %E1, E2, E3
rpmexc = 600;
omegaexc = rpmexc/60*2*pi;
Lelemax=.08; %Longitud máxima de elementos
h = .07;
b = .045;

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
Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
G=E/(2-2*nu);
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
Rrexc=R(isFree);
%% Cuasiestatico
Fcuasi = [0 0 -masapuntual(1,1)*9.81 0 0 0];
Rcuasi = zeros(dof,1);
Rcuasi(n2d6(masa)) = Fcuasi;
Dr = Kr\Rcuasi(isFree);
D = zeros(dof,1);
D(isFree)=Dr;
maxsig=0;
% getTension
% Parametros geometricos cuadrados (Cook pg 50. 2.9-6)
cy=1.5;cz=cy;cT=0.675*h;
%% FORM FUNCTIONS
clear x
xv = [0,1]; % puntos de interpolacion

for e = 1:Nvigas
    Dlocal=D(elemDof(e,:));
    % 1  2  3   4    5    6
    % u  v  w  phi  psi  tita
    u1=Dlocal(1);v1=Dlocal(2);w1=Dlocal(3);ph1=Dlocal(4); p1=Dlocal(5); t1=Dlocal(6);
    u2=Dlocal(7);v2=Dlocal(8);w2=Dlocal(9);ph2=Dlocal(10);p2=Dlocal(11);t2=Dlocal(12);
    n1 = nodos(elementos(e,1),:);
    n2 = nodos(elementos(e,2),:);

    Le=norm(n2-n1);
    ay=12*E*Iy/(Jtors*G*A*Le^2);
    az = 12*E*Iz/(Jtors*G*A*Le^2);
    by=1/(1-ay);
    bz = 1/(1-az);
    sig=0;
    for i=1:length(xv)
        x=xv(i);
        ddv = by*v1*(12*x - 6) - by*v2*(12*x - 6) - (by*t2*(ay - 6*x + 2))/2 + (by*t1*(ay + 6*x - 4))/2;
        ddw = bz*w1*(12*x - 6) - bz*w2*(12*x - 6) - (bz*p2*(az - 6*x + 2))/2 + (bz*p1*(az + 6*x - 4))/2;
        dddv = 3*by*t1 + 3*by*t2 + 12*by*v1 - 12*by*v2;
        dddw = 3*bz*p1 + 3*bz*p2 + 12*bz*w1 - 12*bz*w2;
        Nx = A*E/Le*(u2-u1); 
        Tors = G*K*(ph2-ph1)/Le;
        My = E*Iy*ddw;
        Mz = E*Iz*ddv;
        Vy = E*Iz*dddv;
        Vz = E*Iy*dddw;
        % TENSIONES
        Sx = Nx/A-Mz*h/2/Iz-My*b/2/Iy;
        Txy = Tors*cT/Jtors;
        Ty = cy*Vy/A; Tz = cz*Vz/A;
        
        if Sx>sig
            sig=Sx;
        end
    end
    sig2=0;
    if maxsig<max([sig,sig2]) && max([sig,sig2])<1e9
        maxsig=max([sig,sig2]);
        badel=e;
    end
end

maxsig



resonance

%% 
% for e = 1:Nvigas
%     Dlocal=D(elemDof(e,:));
%     n1 = nodos(elementos(e,1),:);
%     n2 = nodos(elementos(e,2),:);
%     T = Tv(n2-n1,[0 0 1]);
%     Le=norm(n2-n1);
%     [ke, me] = vigastiffness(E,nu,rho,A,Iz,Iy,Jtors,Le);
%     Ke =T*ke*T';
%     
%     flocal=ke*T'*Dlocal;
%     Nx =  flocal(1);
%     Vy =  flocal(2);
%     Vz =  flocal(3);
%     Mx1 = flocal(4);
%     My1 = flocal(5);
%     Mz1 = flocal(6);
%     Mx2 = flocal(10);
%     My2 = flocal(11);
%     Mz2 = flocal(12);
%     Mz=max([Mz1 Mz2]);
%     My = max([My1 My2]);
%     A=b*h;
%     sig = Nx/A+abs(Mz*h/(2*Iz));
%     sig2 = Nx/A+abs(My*b/(2*Iy));
% 
%     if maxsig<max([sig,sig2]) && max([sig,sig2])<1e9
%         maxsig=max([sig,sig2]);
%         badel=e;
%     end
% end

