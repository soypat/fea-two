% DOF  <y> y <theta> 

Ndofpornod = 2;
Nelem = 10;
L=.5; % Regla de 30cm
[nodos3d,elementos,elemDof]=meshViga([0 0 0],[L 0 0],Nelem,Ndofpornod,1);
nodos = nodos3d(:,[1 2]);
% elementos=elem;
%% DOFINITIONS ;)
[Nnod, Ndim] = size(nodos);
dof = Nnod*Ndofpornod;
[Nelem,Ndofporelem]= size(elementos);
% end dofinitions

Le=norm(nodos(1,:)-nodos(2,:));
%viga rectangular. Tomo como ejemplo mi regla de acero
viga=@(E,I,L)(E*I/L^3)*[12   6*L    -12   6*L;
                        6*L  4*L^2  -6*L  2*L^2;
                        -12  -6*L   12    -6*L;
                        6*L  2*L^2  -6*L  4*L^2];
b=.02; %3cm
h=.01;
% h=0.0005; %medio milimetro
I=b*h^3/12;
A=b*h;
E=70e9;
rho=2700;
m = rho*A*Le;
Me= m/420 * [156    22*Le     54      -13*Le
            22*Le    4*Le^2    13*Le    -3*Le^2
            54      13*Le     156     -22*Le
            -13*Le   -3*Le^2   -22*Le   4*Le^2];
Mg = zeros(dof,dof); %lo mismo que zeros(dof)

Ke = viga(E,I,Le);
Kg = zeros(dof,dof);

for e=1:Nelem
    storeTo = elemDof(e,:);
    Kg(storeTo,storeTo) = Kg(storeTo,storeTo)+ Ke;
    Mg(storeTo,storeTo) = Mg(storeTo,storeTo) + Me;
end

% F = v sin wt



isFixed = false(dof,1);
isFixed([2]) = true; % Primer nodo se desplaza libremente en <y>
isFixed([end-1,end]) = true; % Ultimo nodo apoyado simplemente
% isFixed(1:2:end) = true; %fijamos todos los nodos en X para matar modo axial
isFree = ~isFixed;
A = Mg(isFree,isFree)\Kg(isFree,isFree); %F de frecuencia

[Vr, eigVal]=eig(A);
V=zeros(dof);
V(isFree,isFree)=Vr;
Dg =V; % "Desplazamientos globales"
Db = Vr;
%Cantidad de dof despues de reduccion
dofred = size(Vr,2); % Lo mismo que %dofred=sum(isFree)

% PLOTEO ALGO
Dy = Dg(1:2:end,end-1:-2:1); 
n = 1; % Numero de modo.
plot(nodos(1:end,1),Dy(:,n+1))


Phi = zeros(size(Db));
omega = zeros(dofred,1);
omegaray = zeros(dofred,1);
ray = zeros(dofred,1);
for i = 1:dofred
    Dbi=Db(:,i);
    
    aux=Dbi'*Mg(isFree,isFree)*Dbi;
    
    Phi(:,i) = Dbi/sqrt(aux);
    Phii = Phi(:,i);
    
    omegaray(i) = sqrt((Dbi' * Kg(isFree,isFree) *Dbi)/aux);% Cook (11.4-13)
    ray(i)= sqrt((Phii' * Kg(isFree,isFree) *Phii)/(Phii' * Mg(isFree,isFree) *Phii)); %Cook (11.7-1)b
    omega(i) = sqrt(Phii'* Kg(isFree,isFree)*Phii); % Idem
    % Wow son identicas las tres! como es? magia?  <No boludo, es matematica. si M=1 entonces obvio>
end
ESP  = Phi' * Kg(isFree,isFree) *Phi;
omega2 = sqrt(diag(ESP)); % y una cuarta para que tengas

%[omegaray ray omega omega2] %Descomentar para ver que son identicas

%% Damping modal
% La matriz de damping te queda diagonal... Genial! Atentos a la ecuacion
% (11.5-3) del cook
xi1= 0.02; %Damping lo escojo como 0.02 para ambas frecuencias que se yo.
xi2= xi1;
w1 = omega(end);
w2 = omega(end-1);
alfa = 2*w1*w2*(xi1*w2-xi2*w1)/(w2^2-w1^2);
beta = 2*(xi2*w2-xi1*w1)/(w2^2-w1^2);

Cmodal = Phi'*(alfa*Mg(isFree,isFree) +beta*Kg(isFree,isFree))*Phi;
%% Cargas Externas / cargas modales (P y Rmodal [son lo mismo])
Rext = false(dofred,1); % (11.7-5) cook
for i=1:2:size(Mg(isFree,isFree)) % los dof en y los sacudo
    Rext(i,1) = true;
end
Rext= - Mg(isFree,isFree)*Rext;
Rmodal = Phi' *Rext;  %Cada modo se lleva una carga. asi funciona esta wea
% OTRA FORMA:
P = zeros(dofred,1);
for i=1:dofred
    Phii = Phi(:,i);% Ecuacion (11.7-6) cook
    P(i)=Phii' * Rext; % Es lo mismo que Rmodal == P 
end
% El primer modo se lleva Rmodal(1) cantidad de carga 
% Ver descomposicion modal - cargas armonicas
% Lo que estamos haciendo es armar la maquina que va hacer un sine sweep
% al modelo que armamos,

% Si lo hago de la otra forma 1) v= vbar sin wt:

% [Mcc  Mcx]*{D''c} + [Ccc  Ccx]*{D'c} + [Kcc ...]*Dc = 0
% [Mxc  Mxx] {D''x}   [Cxc  Cxx] {D'x}             Dx      % Descomposicion
% conocido/desconocido es una idea loca. hay forma mas facil:
% MxxD''x+CxxD'x + KxxDx = -MxcD''c -CxcD'c -KxcDc
%                        =+Mxcw^2Asen(wt) - CxcAwcos(wt)-KxcAsen(wt)


%% Obtencion de desplazamientos modales
Z = Phi' *Mg(isFree,isFree) ;



