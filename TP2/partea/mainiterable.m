% No-modificables
rojoX = [417 667 1197 1727 2257 2787 3317 3847 4377 4907 5437]'/1000;

ancho = 1.93;
CG = [1.855 ancho/2 1];

%% Mesheador del basamento
meshMe
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
if h>=b
    Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
else
    aux=h;h=b;b=aux;
    Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
    b=h;h=aux;
end
%% Dofinitions
Nnod = size(nodos,1);
n2d6=@(n) [n.*6-5 n.*6-4 n.*6-3 n.*6-2 n.*6-1 n.*6];
n2d3=@(n) [n.*6-5 n.*6-4 n.*6-3];
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

%% ""Rigid Links""
isFixed = false(dof,1);
Nrigid = Nrojo*2;
krigid = 1e8;
rigidElementos =[];
for e = 1:Nrojo %Unimos barras a la masa, para sustituir lo que seria un rigid link
    for i=1:2
    storeTo = [n2d6(rojos(e)) n2d6(masa)];
    n1=nodos(rojos(e,i),:);n2=nodos(masa,:);
    rigidElementos = [rigidElementos;rojos(e,i) masa];
    [ke, ~] = vigastiffness(krigid,nu,0,1,0,0,0,1);
    Ke = vigorotar(ke,n1,n2,[0 0 -1]);
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
    isFixed(storeTo)= isFixed(storeTo) | [0 0 0 0 0 0 0 0 0 1 1 1]'; %Matamos giros en la masa puntual repetidas veces para asegurar su muerte
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
D0 = zeros(dof,1);

Dmasa = [0 0 10e-3 0 0 0]; % se mueve 10mm en y
storeTo = n2d6(masa);
D0(storeTo) = Dmasa;

% Rr = Kr*D0(isFree);
Rr = zeros(size(Kr,1),1);
Rr(n2d3(masa)) = [0 0 4000]';
resonance





