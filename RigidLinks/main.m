%% Rigid Link this One


nodos = [0 0 0
         1 0 0
         2 0 0
         3 0 0
         4 0 0
         1.3 1 1];
masa = 6;
empot= 3;
roller = 5;
rl = [2 4];

elementos = [1 2;2 3;3 4;4 5];
links = [2 6;4 6];

%% Dofinitions
Ndofpornod=6;
[Nelem, Nnodporelem] = size(elementos);
[Nnod, Ndim] = size(nodos);
dof = Ndofpornod*Nnod;
elemDof=node2dof(elementos,Ndofpornod);
n2d = @(n) [n*6-5 n*6-4 n*6-3 n*6-2 n*6-1 n*6];
n2dc = @(n) [n*6-5; n*6-4; n*6-3; n*6-2; n*6-1; n*6];
Nlinks = size(links,1);
%% Vigamaterials
E=200e9;
nu = 0.29;

b=10e-3;
h=5e-3;
A=b*h;
Iz=b*h^3/12;
Iy=h*b^3/12;

Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
rho = 7900;

sz = [0 1 0];
%% Assembler
K = zeros(dof);
for e=1:Nelem
n1 = nodos(elementos(e,1),:);n2 = nodos(elementos(e,2),:);
Le = norm(n2-n1);
storeTo = elemDof(e,:);
[Ke,~] = vigastiffness(E,nu,rho,A,Iz,Iy,Jtors,Le);
[Ke]=vigorotar(Ke,n1,n2,sz);

K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
end

%% Rigids
Clinks = zeros(Nlinks,dof);
aux =[eps eps 1];
for e = 1:Nlinks
    n1 = nodos(links(e,1),:);n2 = nodos(links(e,2),:);
    dof1 = n2d(links(e,1));dof2 = n2d(links(e,2));
    vx = (n2-n1)/norm(n2-n1);    
    Clinks(e,dof1) = [vx(1) vx(2) vx(3) 0 0 0];
    Clinks(e,dof2) = -[vx(1) vx(2) vx(3) 0 0 0];
end
dofL = dof+Nlinks;
Kl = zeros(dofL);
Kl(1:dof,1:dof) = K;
lnk = (1+dof):dofL;
Kl(1:dof,lnk) = Clinks';
Kl(lnk,1:dof) = Clinks;

%% FIXITY
isFixed = false(dofL,1);
isFixed(n2d(empot)) = ones(6,1);
isFixed(n2d(roller)) =  [0 1 1 0 0 0];
isFixed(n2d(masa)) = [0 0 1 1 1 1];
% isFixed(n2d(masa)) = ones(6,1);

isFree = ~isFixed;
% isFixed(lnk) = ones(length(lnk),1);
%% Cargas
R = zeros(dofL,1);
R(n2d(masa)) = [3 2 0 0 0 0]*1e4; % times 10kN
Q=0;
R(dof+1:dofL) = Q;
% Cargas rigid links

%% Mass Locking TEST
isFree(n2d(masa)) = 0;
R(n2d(5)) = [1 1 0 0 0 0]*1e4;
%continue
Krr = Kl(isFree,isFree);
Rrr = R(isFree);
Dr=Krr\Rrr;
D = zeros(dofL,1);
D(isFree) = Dr;

Dnod = [D(1:6:end-Nlinks) D(2:6:end-Nlinks) D(3:6:end-Nlinks)]
