Ge = @(L) 1/L*[-1;1];
kes =@(P,L) P/L*[1 -1;-1 1];
ke  =@(E,A,L) E*A/L*[1 -1;-1 1];
%%Material
E = 200e9;
d = 40e-3; %40mm diametro
A = pi*d^2/4;
%% Dofinitions
nodos = [0 0;1 .4];
elementos = [1 2];
Nelem = size(elementos,1);
Nnod = size(nodos,1);
dof = Nnod*2;
anguloBeta=atand(.4/1)

n2d = @(n) [n*2-1 n*2];
elemDof = [n2d(elementos(:,1)) n2d(elementos(:,2))];
%% Boundary Conditions
isFixed = false(dof,1);
isFixed([1 2 3])= true;
isFree=~isFixed;
%% Cargas;
R = zeros(dof,1);
R(4)=-1;

%% Rigidity
K = zeros(dof,dof);
Ks = zeros(dof,dof);
for e=1:Nelem
    storeTo = elemDof(e,:);
    n1 = nodos(elementos(e,1),:);
    n2 = nodos(elementos(e,2),:);
    Le = norm(n2-n1);
    T = Tbu(atan2d(n2(2)-n1(2),n2(1)-n1(1)));
    Ke = ke(E,A,Le);
    T = Tbu(atan2d(n2(2)-n1(2),n2(1)-n1(1)));
    Ke = T*Ke*T';
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
end
Kr = K(isFree,isFree);
Rr = R(isFree);

Dr = Kr\Rr;
D = zeros(dof,1);
D(isFree)=Dr;

for e = 1:Nelem
    n1 = nodos(elementos(e,1),:);
    n2 = nodos(elementos(e,2),:);
    T = Tbu(atan2d(n2(2)-n1(2),n2(1)-n1(1)));
    Ke = ke(E,A,Le);
    Kerot =T*Ke*T';
    Dlocal=D(elemDof(e,:));
    flocal=Ke*T'*Dlocal;
    %% Acoplo a matriz sigma
    Kes = T*kes(flocal(2),Le)*T';
    Ks(storeTo,storeTo) = Ks(storeTo,storeTo) + Kes;
end
Ksr = Ks(isFree,isFree);
A = Ksr\Kr;
[Vr, lambdacrit] = eig(A);

%Resolución finitos I )



function [T] = Tbu(phi)
%TBU matriz transf. de barra en GRADOS
T=[cosd(phi) 0;
    sind(phi) 0;
    0 cosd(phi);
    0 sin(phi)];
end


