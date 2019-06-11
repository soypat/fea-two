%% Rigid Link this One


nodos = [0 0 0
         1 0 0
         2 0 0
         3 0 0
         4 0 0
         1.3 1 0];
masa = 6;
empot= 3;
roller = 5;
rl = [2 4];

elementos = [1 2;2 3;3 4;4 5];
links = [6 2;6 4];

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

b=50e-3;
h=100e-3;
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

%% Cargas
R = zeros(dof,1);
R(n2d(1)) = [0 1 0 0 0 0]*1e4;

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
Kl=[K Clinks'; Clinks zeros(2)];
% Cargas rigid links
%% FIXITY
isFixed = false(dof,1);
isFixed(n2d(empot))  = ones(6,1);
isFixed(n2d(roller)) =  [0 1 0 0 0 0];
isFixed(n2d(masa))   = [1 1 1 1 1 1];

isFixed(3:6:end) = true;%desp z
isFixed(4:6:end) = true;%giro x
isFixed(5:6:end) = true;%giro y


isFree = ~isFixed;
%% SOLVER
Kr = K(isFree,isFree);
Rr = R(isFree);
Dr=Kr\Rr;
D = zeros(dof,1);
D(isFree) = Dr;
Dnod = [D(1:6:end) D(2:6:end) D(3:6:end)];
posdef = nodos+Dnod;
% Draw_Barra(elementos,nodos,'k')
% Draw_Barra(elementos,posdef,'o')

isFixed(n2d(masa))   = [0 0 1 1 1 1];
isFree = ~isFixed;
dofRelease = [isFree;true(Nlinks,1)];

Klr = Kl(dofRelease,dofRelease);
%% Rigid Links Solve
R(1:dofL) = 0; % Reset Cargas
% R(n2d(1))=[800 0 0 0 0 0]*1e4;% Aumento R en 1
R(n2d(masa)) = [-1 0 0 0 0 0]*1e4; % times 10kN

Q=0;
Rl = [R; Q*ones(Nlinks,1)];
Rlr = Rl(dofRelease);
Dlr=Klr\Rlr;
Dl = zeros(dofL,1);
Dl(dofRelease) = Dlr;
Dlnod = [Dl(1:6:end-Nlinks) Dl(2:6:end-Nlinks) Dl(3:6:end-Nlinks)];
scale=20;
posdefl = nodos+Dlnod*scale;
drawbar3([elementos;links],nodos)
hold on
drawbar3([elementos;links],posdefl)
hold off

% Draw_Barra([elementos;links],posdefl,'k')
% norm(nodos(6,:)-nodos(4,:))
% norm(posdefl(6,:)-posdefl(4,:))

return




