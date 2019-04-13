%%FEA II TP 1 Placas
clear
clc
close all
tic
%% Datos de entrada
%A modificar
Mallado='24X16.xlsx';
CBid='CBa'; %'CBe'  'CBa' %Tipo de cond de borde
tid='t3'; %'t2' 't1' 't3' %Tipo de espesor

E = 210e3;
NU = 0.3;
P=-0.071; %Presion uniforme (MPa) de un doble fondo con calado 7,1m

%Fijos
a=1400; b=1000; % Largo y Ancho
if strcmp(tid,'t1')
    t=a;
elseif strcmp(tid,'t2')
    t=a/10;
elseif strcmp(tid,'t3')
    t=a/100;
end

%% Mallado
A=xlsread(Mallado,'Elementos');
elementos=zeros(size(A,1),size(A,2)-1);
elementos(A(:,1),:)=A(:,2:end);

B=xlsread(Mallado,'Nodos');
nodos=zeros(size(B,1),size(B,2)-1);
nodos(B(:,1),:)=B(:,2:end);


%% Definiciones   
nDofNod = 3;                    % N�mero de grados de libertad por nodo
nNodEle = 4;                    % N�mero de nodos por elemento
nel = size(elementos,1);        % N�mero de elementos
nNod = size(nodos,1);           % N�mero de nodos
nDofTot = nDofNod*nNod;         % N�mero de grados de libertad

%% Condiciones de Borde
bc=false(nNod,nDofNod);
bordes=  nodos(:,1) == 0 | nodos(:,1) == a | nodos(:,2) == 0 | nodos(:,2) == b;
 if strcmp(CBid,'CBe')
     bc(bordes,:)=true;
 elseif strcmp(CBid,'CBa')
     bc(bordes,1)=true;
end
bc=~bc;

 
%%  KIRCHHOFF SYMBOLICA
%% KIRCHHOFF SYMBOLICA
% Asegurarse que todos los elementos generados tienen la misma dimension.
% Mallar en Adina sin biasing o cosas raras
%% Matriz Constitutiva (plane stress)

C = E*t^3/(12*(1 - NU^2))*[ 1.0     NU         0.0
                           NU    1.0         0.0
                          0.0    0.0     (1 - NU)/2 ];
%% Matriz de rigidez
syms x y
func=[1, x, y, x^2, x*y, y^2,x^3,x^2*y,x*y^2,y^3,x^3*y,x*y^3];   

Nodes=nodos(elementos(1,:),:);%Obtencion funciones de forma 
LargoX=(max(Nodes(:,1))-min(Nodes(:,1)))/2;
LargoY=(max(Nodes(:,2))-min(Nodes(:,2)))/2;
Nodes=[-LargoX,-LargoY
        LargoX,-LargoY
        LargoX,LargoY
        -LargoX,LargoY];            
A=sym('A',[12,12]);
for ij=1:4
     A(ij*3-2,:)=subs(func,[x,y],Nodes(ij,:));
     A(ij*3-1,:)=subs(diff(func,x),[x,y],Nodes(ij,:));
     A(ij*3,:)=subs(diff(func,y),[x,y],Nodes(ij,:));
 end
A=double(A);
N=func*A^-1;

Bsym=sym('B',[3,12]);
Bsym(1,:)=diff(diff(N,x),x);
Bsym(2,:)=diff(diff(N,y),y);
Bsym(3,:)=diff(diff(2*N,x),y);

integracionx=[-1,1]*LargoX;
integraciony=[-1,1]*LargoY;
Ke=double(int(int(Bsym'*C*Bsym,x,integracionx),y,integraciony)); %Matriz de rigidez de un elemento
    
K=zeros(nDofTot,nDofTot);
for i=1:nel % Se lo aplica a cada elemento
    eleDofs=node2dof(elementos(i,:),nDofNod);        
    K(eleDofs,eleDofs)=K(eleDofs,eleDofs)+Ke;    
end

%% Carga    
   Q=zeros(12,1);
   Q(1:3:10)=P;
   Re=double(int(int(N'*N*Q,x,integracionx),y,integraciony)); %Fuerza recibida por un elemento
   
Rvec=zeros(nNod*nDofNod,1);
for iele=1:nel % Se o aplica a cada elemento
   eleDofs=node2dof(elementos(iele,:),nDofNod);          
   Rvec(eleDofs)=Rvec(eleDofs)+Re;
end


%% Desplazamientos
Fijo = ~reshape(bc',[],1);
Libre = ~Fijo;
% libre = Libre;
DrK = K(Libre,Libre)\Rvec(Libre);
Dvec = zeros(nDofTot,1);
Dvec(Libre) = Dvec(Libre) + DrK;

DK=reshape(Dvec,[nDofNod,nNod])';

%% Para plot
absisas=min(nodos(:,1)):(nodos(2,1)-nodos(1,1)):max(nodos(:,1));
ordenadas=min(nodos(:,2)):(min(nodos(nodos(:,2) ~= min(nodos(:,2)),2))-nodos(1,2)):max(nodos(:,2));
[X,Y]=meshgrid(absisas,ordenadas);

ZK=reshape(DK(:,1),[size(X,2),size(X,1)]);

%% Recuperacion tensiones

% Descomentar para calcular tensiones con elementos kirchhoff
%GUARDA TARDA MUCHO 750 SEGUNDOS

% uNod=[-LargoX,-LargoY
%         LargoX,-LargoY
%         LargoX,LargoY
%         -LargoX,LargoY];
% stressK=zeros(nNodEle,3,nel);
% for iele = 1:nel
%     nodesEle = nodos(elementos(iele,:),:);
%     for inode = 1:nNodEle
%         
%         A=sym('A',[12,12]);
%         for ij=1:4
%              A(ij*3-2,:)=subs(func,[x,y],Nodes(ij,:));
%              A(ij*3-1,:)=subs(diff(func,x),[x,y],Nodes(ij,:));
%              A(ij*3,:)=subs(diff(func,y),[x,y],Nodes(ij,:));
%          end
%         A=double(A);
%         N=func*A^-1;
% 
%         Bsym=sym('B',[3,12]);
%         Bsym(1,:)=diff(diff(N,x),x);
%         Bsym(2,:)=diff(diff(N,y),y);
%         Bsym(3,:)=diff(diff(2*N,x),y);
% 
%         eleDofs = node2dof(elementos(iele,:),nDofNod);
%         stressK(inode,:,iele)=C*double(subs(subs(Bsym,x,uNod(inode,1)),y,uNod(inode,2)))*Dvec(eleDofs);
%     end
% end





%% MINDLIN ISOPARAMETRICO
%% MINDLIN ISOPARAMETRICO
%% Matrices Constitutivas
% De flexion
Db=E*t^3/(12*(1-NU^2));
Cb = Db*[ 1.0     NU         0.0
                           NU    1.0         0.0
                          0.0    0.0     (1 - NU)/2 ];
                      
% De corte
Cs=(5/6)*t*E/(2*(1+NU))*eye(2);

%% Matriz de rigidez
K=zeros(nDofTot);

% Matriz Bending
rsInt = 2*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

for iele=1:nel
    Ke=zeros(nDofNod*nNodEle);
    nodesEle=nodos(elementos(iele,:),:);
    for ipg=1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],'Q4');
        Nf=N;
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4');
        dNf=dN;
        % Derivadas de x,y, respecto de ksi, eta
        nodf=nodesEle;
        jac = dN*nodesEle;
        Djac=det(jac);
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        %B Bending Matrix Bb
        Bb=zeros(3,12);
        Bb(1,2:3:11)=dNxy(1,:);
        Bb(2,3:3:12)=dNxy(2,:);
        Bb(3,2:3:11)=dNxy(2,:);
        Bb(3,3:3:12)=dNxy(1,:);
        Bbf=Bb;
        Ke=Ke+Bb'*Cb*Bb*wpg(ipg)*Djac;
        Kef=Ke;
    end
    eleDofs=node2dof(elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs)=K(eleDofs,eleDofs)+Ke;
end

% Matriz shear
rsInt = 1*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);

for iele=1:nel
    Ke=zeros(nDofNod*nNodEle);
    nodesEle=nodos(elementos(iele,:),:);
    for ipg=1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Funciones de forma respecto de ksi, eta
        N = shapefuns([ksi eta],'Q4');
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        Djac=det(jac);
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;
        
        %B shear Bs
        Bs=zeros(2,12);
        Bs(1,1:3:10)=-dNxy(1,:);
        Bs(1,2:3:11)=N;
        Bs(2,1:3:10)=-dNxy(2,:);
        Bs(2,3:3:12)=N;
        
        Ke=Ke+Bs'*Cs*Bs*wpg(ipg)*Djac;
    end
    eleDofs=node2dof(elementos(iele,:),nDofNod);
    K(eleDofs,eleDofs)=K(eleDofs,eleDofs)+Ke;
end
Kf = K;

%% Carga

F=zeros(nDofTot,1);
rsInt = 1*ones(1,2);
[wpg, upg, npg] = gauss(rsInt);
for iele=1:nel
    nodesEle=nodos(elementos(iele,:),:);
    
    Area=0;
    for ipg=1:npg
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder([ksi eta],'Q4');
        % Derivadas de x,y, respecto de ksi, eta
        jac = dN*nodesEle;
        Djac=det(jac);
          
        Area=Area+ wpg(ipg)*Djac;  %Area del elemento            
    end
    carga=P*Area;
    eleDofs=node2dof(elementos(iele,:),nDofNod);
    F(eleDofs(1:3:10))=F(eleDofs(1:3:10))+carga/4*ones(4,1);
end
Rf=F;
%% Desplazamientos
Fijo = ~reshape(bc',[],1);
Libre = ~Fijo;
libre=Libre;
DrM = K(Libre,Libre)\F(Libre);


Dvec = zeros(nDofTot,1);
Dvec(Libre) = Dvec(Libre) + DrM;

DM=reshape(Dvec,[nDofNod,nNod])';
ZM=reshape(DM(:,1),[size(X,2),size(X,1)]);
% Grafico PATO
D=Dvec;
Df=D;
Dzf = zeros(25,17); % Matriz superficie
xv =[];
yv = [];
Nnod =25*17;
for n=1:Nnod
    xv = [xv nodos(n,1)];
    yv = [yv nodos(n,2)];
    Dzf(n) = D(n*3-2);
end
Dzf=reshape(Dzf,[],1)';
scatter3(xv,yv,Dzf)

%% Recuperacion de tensiones en los nodos
uNod = [-1 -1
             1 -1
             1  1
            -1  1];
 
stressb = zeros(nNodEle,3,nel);
stressS = zeros(nNodEle,2,nel);
for iele = 1:nel
    nodesEle = nodos(elementos(iele,:),:);
    for inode = 1:nNodEle
        % Punto de Gauss
        ksi = uNod(inode,1);
        eta = uNod(inode,2);
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN  = shapefunsder([ksi eta],'Q4');
        % Derivadas de x,y, respecto d % '2'; % e ksi, eta
        jac  = dN*nodesEle;
        % Derivadas de las funciones de forma respecto de x,y.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        Bb=zeros(3,12);
        Bb(1,2:3:11)=dNxy(1,:);
        Bb(2,3:3:12)=dNxy(2,:);
        Bb(3,2:3:11)=dNxy(2,:);
        Bb(3,3:3:12)=dNxy(1,:);
        
        Bs=zeros(2,12);
        Bs(1,1:3:10)=-dNxy(1,:);
        Bs(1,2:3:11)=N;
        Bs(2,1:3:10)=-dNxy(2,:);
        Bs(2,3:3:12)=N;
        
        eleDofs = node2dof(elementos(iele,:),nDofNod);
        
        stressb(inode,:,iele) = Cb*(Bb*Dvec(eleDofs));
        stressS(inode,:,iele) = Cs*(Bs*Dvec(eleDofs));
    end
end



%% Analitica
%% Analitica
if strcmp(CBid,'CBa')
    W_an=  @(x,y) 0;
    for n=1:2:21
        for m=1:2:21
            W_an= @(x,y) W_an(x,y) + 16*P/(pi()^6*Db).*(sin(m.*pi().*x./a).*sin(n.*pi().*y./b))./(m*n*((m/a)^2+(n/b)^2)^2);
        end
    end 
    z_an=W_an(nodos(:,1),nodos(:,2));
    z_anbc=z_an;
    z_an(~bc(:,1))=0;
    z_anbc(~bc(:,1))=1;
    
    Zex=reshape(z_an,[size(X,2),size(X,1)]);
    Zexbc=reshape(z_anbc,[size(X,2),size(X,1)]);

elseif strcmp(CBid,'CBe')
%RESPUESTA ANALITICA DE EMPOTRADA
end




%% Ploteo
figure (1)
meshplot(elementos,nodos,'b', 1)

figure(2)
mesh(X,Y,ZK')
title('Deflexion Kirchhoff')

figure(3)
mesh(X,Y,ZM')
title('Deflexion Mindlin')

figure(4)
mesh(X,Y,abs(Zex'-ZK')./abs(Zexbc')*100)
title('Error % Deflexion Kirchhoff')

figure(5)
mesh(X,Y,abs(Zex'-ZM')./abs(Zexbc')*100)
title('Error % Deflexion Mindlin')

figure(6)
bandplot(elementos,nodos,squeeze(stressb(:,1,:))',[],'k') %bandplot(elementos,nodos,squeeze(stressS(:,1,:))',[],'k')
title('Tensiones Mindlin')
Sbf=stressb;
save('flo','Dzf','Df','Kf','Rf','libre','Sbf');
% figure(7)
% bandplot(elementos,nodos,squeeze(stressK(:,1,:))',[],'k')
% title('Tensiones Kirchhoff')


toc
