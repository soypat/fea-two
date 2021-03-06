%Placa Kirchoff
funcForma
N=shapefuns;

%% Problema - nodos/elementos

divisionesx = 4;
divisionesy = 4;
a=1;b=1; % Tama�o del problema
dx = a/(divisionesx-1);
dy = b/(divisionesy-1);

[nodos, elementos] = rectmesh(a,b,divisionesx,divisionesy);
%% DOFINITIONS
[Nelem, Nnodporelem]= size(elementos);  
[Nnod, Ndim] = size(nodos); % Numero de nodos, Numero de dimensiones del problema (es 2-D)

Ndofpornod = 3; %Para placa Kirchoff
dof = Nnod*Ndofpornod;
DOF = (1:dof)'; %vector columna

n2d = @(nodo) [nodo*3-2, nodo*3-1, nodo*3]; % Funci�n Node a DOF. Obtiene indices de dof de un nodo. Si hay mas/menos de 3 dof por nodo entonces cambia 
%% Bonus: Armo matriz elemDof para ensamblar Matriz rigidez rapido
elemDof = zeros(Nelem,Ndofpornod*Nnodporelem);
for e = 1:Nelem
   for n = 1:Nnodporelem
       elemDof(e,n2d(n)) = n2d(elementos(e,n));
   end
end

%% Condiciones de Borde y Cargas
isFixed = false(dof,1);
% simplementeApoyadoNodos = [1 2 3 4 5 10 15 20 25 24 23 22 21 16 11 6]';
for n =1:Nnod
    if nodos(n,1)==0 || nodos(n,2)==0 || nodos(n,2)==b || nodos(n,1)==a
        isFixed(n2d(n))=true;
    end
end
isFree = ~isFixed;

syms x y real
% q = sin(pi*x)*sin(pi*y); % Caso particular para el problema: LA carga == Los desplazamientos. Ver placas-> Ugural.

%% Propiedades del material
E = 1.000;
NU = 0.300;
t = 1.0;
F = E*t^3/(12*(1 - NU^2)); % Rigidez a la flexion.
C = F*[ 1.0     NU         0.0
                    NU    1.0         0.0
                   0.0    0.0     (1 - NU)/2 ];

%% Obtencion Matriz Rigidez por integracion simbolica No-Isoparametrica
Kg = sparse(dof,dof);
Ke = double(int(int(B'*C*B,x,0,a/(divisionesx-1)),y,0,b/(divisionesy-1))); %Los elementos son 1x1. Integro asi
for e = 1:Nelem
    storeTo = false(dof,1);
    storeTo(elemDof(e,:)) = true;
%     storeTo = elemDof(e,:);
    Kg(storeTo,storeTo) = Kg(storeTo,storeTo) + Ke;
%     eval(sum(diag(Ke)))
end

%% Obtencion vector columna de cargas R
%Para el caso que quiero b=2, a=2.
p0=1;
q = p0;%p0*sin(x*pi/a)*sin(y*pi/b);%*sin(pi*x/2)*sin(pi*y/2); % Caso particular para el problema: LA carga == Los desplazamientos. Ver placas-> Ugural.
R = zeros(dof,1);

e=1; %Elementos mismo tama�o
    x1 = min(nodos(elementos(e,:),1));
    x2 = max(nodos(elementos(e,:),1));
    y1 = min(nodos(elementos(e,:),2));
    y2 = max(nodos(elementos(e,:),2));
A = (x2-x1)*(y2-y1);
% for e = 1:Nelem %FORMA Rapida
%     storeTo = false(dof,1);
%     storeTo(elemDof(e,:))=true;
%     x =(max(nodos(elementos(e,:),1))- min(nodos(elementos(e,:),1)))/2+min(nodos(elementos(e,:),1));
%     y = (max(nodos(elementos(e,:),2))-min(nodos(elementos(e,:),2)))/2+min(nodos(elementos(e,:),2));
%     
%     
%     Q = eval(N.'*p0*sin(x*pi/a)*sin(y*pi/b)*A);
%     R(storeTo)=R(storeTo)+Q;
% end

Nint = double(int(int(N.',x,0,dx),y,0,dy));
for e = 1:Nelem %FORMA EXACTA MUY LENTA
    storeTo = false(dof,1);
    storeTo(elemDof(e,:))=true;
%     x1 = min(nodos(elementos(e,:),1));
%     x2 = max(nodos(elementos(e,:),1));
%     y1 = min(nodos(elementos(e,:),2));
%     y2 = max(nodos(elementos(e,:),2));
%     Q = int(int(N.'*q,x,x1,x2),y,y1,y2);
    Q = Nint*p0;
    R(storeTo)=R(storeTo)+Q;
end

%% Inversion (C�lculo de desplazamientos)
Dr = Kg(isFree,isFree)\R(isFree);


D=zeros(dof,1);
D(isFree) = Dr;

Dz = zeros(divisionesx,divisionesy); % Matriz superficie
xv =[];
yv = [];
for n=1:Nnod
    xv = [xv nodos(n,1)];
    yv = [yv nodos(n,2)];
    Dz(n) = D(n*3-2);
end
Dz=reshape(Dz,[],1)';
scatter3(xv,yv,Dz)
xv =[]; yv=[];
qn = [];
for n = 1:Nnod
    x=nodos(n,1); y = nodos(n,2);
    xv = [xv nodos(n,1)]; yv = [yv nodos(n,2)];
    qn = [qn p0*sin(x*pi/a)*sin(y*pi/b)];
    
end
qn=double(qn);
figure
scatter3(xv,yv,qn);


    
% %% Solucion Exacta
% Niter = 50;
% w_max=0;
% for m = 1:2:Niter
%    for n = 1:2:Niter %Ugural 4ta ed. Ecuacion (13.22)
%        w_max=w_max + ((-1).^((m+n)/2-1))/((m*n* ( (m/a)^2+(n/b)^2 )^2)); 
%    end
% end
% 
% w_max = 16*p0/(pi^6)./F*w_max;
% fprintf('Para un carga uniforme el desplazamiento maximo es\nw_max=%f\n',w_max)
% w=D(13*3-2);
% fprintf('A vos te dio w_max=%f',w)
%% Secci�n debug
% debug=false;
% if debug==false
% %     return
% end
% fprintf('La matriz rigidez reducida:\n')
% disp(Kg(isFree,isFree));
% % fprintf('Las cargas reducidas:\n')
% % disp(R(isFree));
% fprintf('Los desplazamientos:\n')
% disp(D(n2d(13))) %Mi nodo central es el nodo de la mala suerte. Mala eleccion por mi parte
% 
% elPloto(nodos,elementos)
