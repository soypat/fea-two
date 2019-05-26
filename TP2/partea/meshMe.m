%% Código
soporteX = [0 longE(1)+.097 longE(2)+longE(1)+.097 longE(1)+longE(2)+longE(3)+.097 LargoMotor]';
azulX = [0.097 AB(1)+.097 sum(AB)+.097 5.567];
Nazul = length(azulX);
Nrojo = length(rojoX);
Nsop = length(soporteX);
elementos =[];

%% Mesheo primer soporte
nodos = [0 0 0;0 ancho 0];
[nodos, elementos] = mesh1DL(nodos,1,2,Lelemax);
lastTop = 2;
lastBot = 1;
%% Voy a guardar a memoria cuales son mis nodos
rojos = zeros(Nrojo,2);
azules = zeros(Nazul,2);
soportes = zeros(Nsop,2);

krojo= 1;
ksop = 2;
kazul=1;
x=0.00;

keepGoing =true;

while keepGoing
    if ksop>Nsop
        keepGoing=false;
        break
    end
    caseNode = nextNode(x,soporteX,rojoX,azulX);
    Nnod = size(nodos,1);
    switch caseNode
        case 1 %Caso soporte
            xnext=soporteX(ksop);
            nodos = [nodos;xnext 0 0;xnext ancho 0];
            soportes(ksop,:) = [Nnod+1 Nnod+2];
            [nodos, elem1] = mesh1DL(nodos,Nnod+1,Nnod+2,Lelemax); %Mesheo el soporte de abajo para arriba
            [nodos, elem2] = mesh1DL(nodos,lastBot,Nnod+1,Lelemax); %Mesheo viga longit. bottom
            [nodos, elem3] = mesh1DL(nodos,lastTop,Nnod+2,Lelemax); % Mesheo """ top
            elementos = [elementos;elem1;elem2;elem3];
            ksop = ksop+1;
        case 2 %Caso rojo - Motor a BASAMENTO - Rigid linkers. Fijos
            xnext = rojoX(krojo);
            nodos = [nodos;xnext 0 0;xnext ancho 0];
            rojos(krojo,:) = [Nnod+1 Nnod+2];
            [nodos, elem1] = mesh1DL(nodos,lastBot,Nnod+1,Lelemax); %Mesheo viga longit. bottom
            [nodos, elem2] = mesh1DL(nodos,lastTop,Nnod+2,Lelemax); % Mesheo """ top
            elementos = [elementos;elem1;elem2];
            krojo=krojo+1;
        case 3
            xnext = azulX(kazul);
            nodos = [nodos;xnext 0 0;xnext ancho 0];
            azules(kazul,:) = [Nnod+1 Nnod+2];
            [nodos, elem1] = mesh1DL(nodos,lastBot,Nnod+1,Lelemax); %Mesheo viga longit. bottom
            [nodos, elem2] = mesh1DL(nodos,lastTop,Nnod+2,Lelemax); % Mesheo """ top
            elementos = [elementos;elem1;elem2];
            kazul=kazul+1;
    end
    lastBot = Nnod+1;
    lastTop = Nnod+2;
    x=xnext;
end

