clear
close all
clc
tic
%% Datos

c = 0.02;
E = 0.7e11; % Pa
rho = 2700; % Kg/m^3
A = c^2; I = c^4/12;
L=0.5; % m

%% Discretizacion

nNod=9; % Numero de nodos
nDofNod = 2; % Numero de dof por nodos

nodos = [linspace(0,L,nNod)' zeros(nNod,1)];
elementos = [(1:nNod-1)' (2:nNod)'];

nEle = size(elementos,1);
nDofTot = nNod*nDofNod;

dof=reshape(1:nDofNod*nNod,[nDofNod,nNod])';

%% Condiciones de borde
bc = false(nNod,nDofNod);
bc(1,:) = true;

isFixed = reshape(bc',[],1);
isFree = ~isFixed;

%% Matriz de Rigidez y de Masa
K = zeros(nDofTot);
M = zeros(nDofTot);

for iele=1:nEle
    v=nodos(elementos(iele,2),:)-nodos(elementos(iele,1),:);
    Lele=norm(v);
    v=v/Lele;
    m=A*Lele*rho;
    
    Y1 = 12*E*I/Lele^3; Y2 = 6*E*I/Lele^2; Y3 = 4*E*I/Lele; Y4 = 2*E*I/Lele;
    
    kele=[Y1 Y2 -Y1 Y2;Y2 Y3 -Y2 Y4;-Y1 -Y2 Y1 -Y2;Y2 Y4 -Y2 Y3];
    
    mele = m/420*[156 22*Lele 54 -13*Lele; 22*Lele 4*Lele^2 13*Lele -3*Lele^2;...
        54 13*Lele 156 -22*Lele; -13*Lele -3*Lele^2 -22*Lele 4*Lele^2];    
    
    ubic=[dof(elementos(iele,1),:) dof(elementos(iele,2),:)];
    
    K(ubic,ubic)=K(ubic,ubic)+kele;
    M(ubic,ubic)=M(ubic,ubic)+mele;
end

%% Obtencion de desplazamiento y frecuencia
[eigVect,eigVal] = eig(M(isFree,isFree)\K(isFree,isFree));
eigVect = flip(eigVect,2); eigVal = flip(flip(eigVal,1),2);

W=sqrt(diag(eigVal));

%% Plot

% Nmodo=2; % Nro de modo    
% desp=zeros(1,nDofTot);
% desp(isFree)=eigVect(:,Nmodo);
% Nro_ptos_ele=15;
% desp_tot =0; X=0;
% 
% for iele=1:nEle
%     
%         v=nodos(elementos(iele,2),:)-nodos(elementos(iele,1),:);
%         Lele=norm(v);
%     
%         x0 = nodos(elementos(iele,1),1);
%         xL = nodos(elementos(iele,2),1);
% 
%         x = linspace(x0,xL,Nro_ptos_ele);
%         X=[X,x(2:end)];
% 
%         N = [1 - 3*(x-x0).^2/Lele^2 + 2*(x-x0).^3/Lele^3;
%             (x-x0) - 2*(x-x0).^2/Lele + (x-x0).^3/Lele^2;
%             3*(x-x0).^2/Lele^2 - 2*(x-x0).^3/Lele^3;
%             -(x-x0).^2/Lele + (x-x0).^3/Lele^2]';
%         
%         pos = [dof(elementos(iele,1),:) dof(elementos(iele,2),:)];
% 
%         desp_ele = N*desp(pos)';
%         desp_tot = [desp_tot; desp_ele(2:end)];
% end
% 
% plot(X,desp_tot,'b')
% hold on
% plot(nodos(:,1),desp(1:2:end-1),'or')    
% grid on
%% Descomposicion Modal

D=zeros(nDofTot,nDofTot-sum(isFixed));
D(isFree,:)=eigVect;
Phi=D;

for i_D=1:size(D,2) %Normalizado en base a M
    
    modulo=(D(:,i_D)'*M*D(:,i_D));  
    Phi(:,i_D)=D(:,i_D)./sqrt(modulo);
    
end

%% Calculo amplitud en funcion de frecuencias
Rango_Omega=0:20:10000;
Rango_Ksi=0.05:0.05:0.3;

acceleracion=zeros(nDofTot,1); 
acceleracion(1)=1;
R_ext=-M*acceleracion;

R_Phi=Phi'*R_ext;
R_Phi=abs(R_Phi); % Esto no es robado???

n_Omega=1;
Amplitudes=zeros(length(Rango_Omega),length(Rango_Ksi));
for Omega=Rango_Omega
    n_Ksi=1;
    for Ksi=Rango_Ksi
        z=zeros(length(W),1);
        for i=1:length(W)
            Beta=Omega/W(i);
            z(i)=(R_Phi(i)/(W(i)^2))/sqrt((1-Beta^2)^2+(2*Ksi*Beta)^2);            
        end
        Amplitudes(n_Omega,n_Ksi)=sum(z);
        n_Ksi=n_Ksi+1;
    end     
    n_Omega=n_Omega+1;
end

%Amplitudes=Amplitudes*10^9;

semilogy(Rango_Omega,Amplitudes(:,1),'r'); 
hold on
semilogy(Rango_Omega,Amplitudes(:,2),'g');
hold on
semilogy(Rango_Omega,Amplitudes(:,3),'b'); 
hold on
semilogy(Rango_Omega,Amplitudes(:,4),'k'); 
hold on
semilogy(Rango_Omega,Amplitudes(:,5),'y');
hold on
semilogy(Rango_Omega,Amplitudes(:,6),'m'); 
grid on

% plot(Rango_Omega,Amplitudes(:,1),'ro'); 
% hold on
% plot(Rango_Omega,Amplitudes(:,2),'go');
% hold on
% plot(Rango_Omega,Amplitudes(:,3),'bo'); 
toc