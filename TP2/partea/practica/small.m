
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

for e = 1:Nvigas
    storeTo = elemDof(e,:);
    n1=nodos(elementos(e,1),:);n2=nodos(elementos(e,2),:);
    L=norm(n2-n1);
    ke = vigastiffness(E,nu,A,Iz,Iy,Jxy,L);
    Ke = vigorotar(ke,n1,n2,[0 0 -1]);
    K(storeTo,storeTo) = K(storeTo,storeTo) + Ke;
end


