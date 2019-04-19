Ndofpornod = 2;
Nelem = 5;
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
h=b;
% h=0.0005; %medio milimetro
I=b*h^3/12;
A=b*h;
E=.7e11;
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

isFixed = false(dof,1);
isFixed([1 2]) = true; %fijo los primeros 3 dof
% isFixed(1:2:end) = true; %fijamos todos los nodos en X para matar modo axial
isFree = ~isFixed;
A = Mg(isFree,isFree)\Kg(isFree,isFree); %F de frecuencia

[Vr, eigVal]=eig(A);
V=zeros(dof);
V(isFree,isFree)=Vr;
Vr=Vr(:,end:-1:1); %X_i
w2 = diag(eigVal);
w2= w2(end:-1:1);% w_i^2
w=w2.^.5;     %w_i
hz = w/(2*pi); %frecuencia

%% Plot mode 1
mode = 1;
elediv = 10; %Indica que tan refinado es el analisis de desplazamiento
NmicroElem = Nelem*(elediv-1);
xtot = linspace(0,L,NmicroElem+1);
desp = [0 0 Vr(:,mode)']; %Que es esto????
Detot=0;

for e = 1:Nelem
    elemNod = nodos(elementos(e,:),:);
    Le = norm(elemNod(end,:)-elemNod(1,:));
    x0 = elemNod(1,1);
    xe = linspace(x0,elemNod(end,1),elediv);
    
    N = [1 - 3*(xe-x0).^2/Le^2 + 2*(xe-x0).^3/Le^3
            (xe-x0) - 2*(xe-x0).^2/Le + (xe-x0).^3/Le^2
            3*(xe-x0).^2/Le^2 - 2*(xe-x0).^3/Le^3
            -(xe-x0).^2/Le + (xe-x0).^3/Le^2]';
    De = N*desp(elemDof(e,:))';
    Detot = [Detot;De(2:end)];%desplazamiento del elemento
end
plot(xtot,Detot)
% rayleigh = zeros(size(eigVal,1),1);
% for j=1:length(w)
%    xj = Vr(:,j);
%    xi = xj;
%    kx=xi'*Kg(isFree,isFree)*xj; % si i==j no dan cero 
%    mx=xi'*Mg(isFree,isFree)*xj;  
%    rayleigh(j) = kx/mx;
% end

% Observar que w2==rayleigh (muy similares)
%% Matriz Espectral
W2=diag(rayleigh);

%% Ensamble estricto de D_raya
%Pasamos a notación desplazamientos (sin angulo/segundo dof)
D = V(1:2:end,:);
isFreeD=isFree(1:2:end);

for i = 1:4
    Di=D(:,i);
    Dj=Di;
%     Di'*Mg()
end
% figure
% for p=1:2:5
%     plot(nodos(:,1),D(:,p))
%     hold on
% end
% figure
% for p=2:2:6
%     plot(nodos(:,1),D(:,p))
%     hold on
% end

return
nodos = 12*[-5,0;-5,10;5,10;5,0];%in
E=30e6;%psi
A=10;%in^2
I13=200;%in^4
I2=100;
L=10*12;
% viga=@(E,I,L)(E*I/L^3)*[12 6*L -12 6*L;
%     6*L 4*L^2 -6*L 2*L^2;
%     -12 -6*L 12 -6*L;
%     6*L 2*L^2 -6*L 4*L^2];
vigota=@(E,A,I,L) [A*E/L 0 0 -A*E/L 0 0;
    0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
    0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
    -A*E/L 0 0 A*E/L 0 0;
    0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
    0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L];

Tvu=@(phi) [cosd(phi) sind(phi) 0 0 0 0;
    -sind(phi) cosd(phi) 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 cosd(phi) sind(phi) 0;
    0 0 0 -sind(phi) cosd(phi) 0;
    0 0 0 0 0 1];
Ne=3;
N=Ne+1;
ndof=3; %dof unitarios
Ndof=N*ndof;
elementos=zeros(Ne,6);
for i = 1:Ne
    for j = 1:6
        elementos(i,j)=(i-1)*ndof+j;
    end
end



kG=zeros(Ndof);

T90=Tvu(90);
T0=Tvu(0);
Tm90=Tvu(-90);

k13=vigota(E,A,I13,L);
k=vigota(E,A,I2,L);
k1=T90'*k13*T90;
k2=T0'*k*T0;
k3=Tm90'*k13*Tm90;

kG(elementos(1,:),elementos(1,:))=kG(elementos(1,:),elementos(1,:))+k1;
kG(elementos(2,:),elementos(2,:))=kG(elementos(2,:),elementos(2,:))+k2;
kG(elementos(3,:),elementos(3,:))=kG(elementos(3,:),elementos(3,:))+k3;

CB=false(Ndof,1);
CB(1)=true;
CB(2)=true;
CB(3)=true;
CB(end-2)=true;
CB(end-1)=true;
CB(end)=true;

K=kG(~CB,~CB);



R=zeros(Ndof,1);
R(4)=10e3;
R(9)=5e3;
F=R(~CB);
U=K\F;%
D=zeros(ndof*ndof);
D(~CB) = U;