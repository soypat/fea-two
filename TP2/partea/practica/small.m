
rojoNodos = [1 0 0;1.2 0 0];
azulNodos = [0 0 0;3 0 0];

Nrojo = size(rojoNodos,1);
Nazul = size(azulNodos,1);

rojos = zeros(Nrojo,1); % apunta al nodo correspondiente
azules = zeros(Nazul,1);


nodos=[];
elementos =[];
% nodos = [0 0 0];
Nnod = 0;

keepGoing =true;
lookingForSecondNode = false;
elNodo = 1;
nodo = azulNodos(1,:);
nodos = nodo;
Nele =1;
krojo= 1;
kazul=2;
azules(1)=1;

while keepGoing %CADA ITERACION DEBERIA CREAR UN TRAMO
    if krojo>Nrojo
        nodos = [nodos; azulNodos(kazul,:)]; %Mesheamos el nodo azul
        azules(kazul) = size(nodos,1);
        % MESHEAMOS TRAMO
        [nodos,elem] = mesh1D(nodos,elNodo,azules(kazul),Nele);
        elementos= [elementos;elem];
        elNodo = size(nodos,1);
        break
    end
    %% Encontramos el segundo nodo del tramo y mesheamos
    if rojoNodos(krojo,1)<azulNodos(kazul,1)
        nodos = [nodos; rojoNodos(krojo,:)]; %Mesheamos el nodo rojo 
        rojos(krojo) = size(nodos,1);
        % MESHEAMOS TRAMO
        [nodos,elem] = mesh1D(nodos,elNodo,rojos(krojo),Nele);
        elNodo = rojos(krojo);
        krojo=krojo+1;
    elseif rojoNodos(krojo,1)>azulNodos(kazul,1)
        nodos = [nodos; azulNodos(kazul,:)]; %Mesheamos el nodo azul
        azules(kazul) = size(nodos,1);
        % MESHEAMOS TRAMO
        [nodos,elem] = mesh1D(nodos,elNodo,azules(kazul),Nele);
        kazul=kazul+1;
         elNodo = azules(kazul);
    end
    elementos= [elementos;elem];
end

% nodoMasa=[1.1 2 0];
% 
% nodos=[nodos;nodoMasa];

masa=size(nodos,1);
Nvigas = size(elementos,1);
%% Vigamaterials
E = 200e9;
rho = 7900;
nu=0.3;
h=.02;
b= .02;
A=b*h;
Iz = b*h^3/12;
Iy = h*b^3/12;
Jxy = h*b*h*b/64;





%% Dofinitions

Nnod = size(nodos,1);
n2d6=@(n) [n.*6-5 n.*6-4 n.*6-3 n.*6-2 n.*6-1 n.*6];
elemDof = [n2d6(elementos(:,1)) n2d6(elementos(:,2))];
dof = Nnod*6;
K = zeros(dof,dof);
M = zeros(dof,dof);

for e = 1:Nvigas
    storeTo = elemDof(e,:);
    n1=nodos(elementos(e,1),:);n2=nodos(elementos(e,2),:);
    Le=norm(n2-n1);
  
    [ke, me] = vigastiffness(E,nu,rho,A,Iz,Iy,Jxy,Le);
    Ke = vigorotar(ke,n1,n2,[0 0 1]);
    Me = vigorotar(me,n1,n2,[0 0 1]);
    
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
    M(storeTo,storeTo) = M(storeTo,storeTo) + Me;
end
%% Comenzamos a definir las condiciones de borde así matamos giros en barras
isFixed = false(dof,1);
Nrigid = Nrojo;
for e = 1:Nrigid %Unimos barras a la masa, para sustituir lo que seria un rigid link
    storeTo = [n2d6(rojos(e)) n2d6(masa)];
    n1=nodos(rojos(e),:);n2=nodos(masa,:);
    Le=norm(n2-n1);
    [ke, ~] = vigastiffness(E,nu,0,A*10,0,0,0,Le);
    Ke = vigorotar(ke,n1,n2,[0 0 -1]);
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
    isFixed(storeTo)= isFixed(storeTo) | [0 0 0 0 0 0 0 0 0 1 1 1]'; %Matamos giros en la masa puntual repetidas veces para asegurar
end
%% Acoplo masa puntual
masapuntual = 1e6*blkdiag(.055, .055, .055, 0.49, 0.28455,0.28455);
storeTo = n2d6(masa);
M(storeTo,storeTo)=M(storeTo,storeTo) + masapuntual;
%% Mato z sobre masa puntual
storeTo = n2d6(masa);
isFixed(storeTo(3))=true;
%% Mato giro en x sobre viga
dofViga = reshape(n2d6([rojos;azules])',[],1);
isFixed(dofViga(4:6:end)) = true;

%% Creo las condiciones de bordes en forma de bulones
Eb=200e9/1e3;
d =20e-3; %20mm
Ab= d^2/4*pi;
Lb =70e-3; %70mm
Lbtransversal = d;
kb = Eb*Ab/Lb;
kbt = Eb*Ab/Lbtransversal;
Kg = abulonar(azules,K,[kb kbt kbt]); %Hice esta funcion pesimamente. devuelve tamaño inecesario
%% Aplico condiciones de borde! En este caso es impedir que giren las barras conectadas a la masa puntual
K = Kg(1:dof,1:dof);
isFree=~isFixed;
Kr = K(isFree,isFree);
Mr = M(isFree,isFree);
digitosPerdidos = get_cond(Kr);
if digitosPerdidos >10
    error('Condicionamiento malo!')
end


resonance
