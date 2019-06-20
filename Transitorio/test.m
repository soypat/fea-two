
Ndofpornod = 2;
Nelem = 1;
L=.5; % 0.5 metros % TODO UNIDADES ESTANDAR
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
b=20e-3; % 3cm
h=10e-3; % 1cm
% h=0.0005; %medio milimetro
I=b*h^3/12;
Amp=b*h;
E=70e9; %Pa
rho=2700;
m = rho*Amp*Le;
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
isFixed([1 2]) = true; % Primer nodo empotrado
isFree = ~isFixed;
K = Kg(isFree,isFree);
M = Mg(isFree,isFree);

dofred = sum(isFree);
A = Mg(isFree,isFree)\Kg(isFree,isFree); %F de frecuencia
[Vr, eigVal] = eig(A);
Phi = zeros(dofred,dofred);
omega = zeros(dofred,1);
for i = 1:dofred
   Dbi=Vr(:,i);
   Phi(:,i) = Dbi/sqrt(Dbi'*M*Dbi);
   Phii = Phi(:,i); 
   omega(i) = sqrt(Phii'* Kg(isFree,isFree)*Phii); % Idem
end
alpha = 0.02;
beta = 0.0002%.1;
C = alpha*M+beta*K;
% cmod=Phi'*(alpha*Mr+beta*Kr)*Phi

dt = 1/max(omega)/5000; %segundos CONDICION: dt <= 2/max(omega)
dtdiv2 = dt/2;
dtdiv6 = dt/6;
Nt = 2e7; % 10 segundos simulados
% Minv = inv(Mg(isFree,isFree));
% Kinv = inv(Kg(isFree,isFree));
%% Cargas
Rext = zeros(dof,1);
Rext2 = zeros(dof,1);
q=1;
freqCarga = 1000; %Hz

V = zeros(dofred,1); % Condicion inicial: velocidad = 0
V(end-1) = 5;
D = zeros(dofred,1); % desplazamientos = 0
D_t = zeros(dofred,Nt);
for i = 0:Nt
    t = dt*i; %tiempo
    if true %t<=dt*100
        Rext(end)=0;%q*sin(freqCarga*t*2*pi+eps);%20kN
        Rext2(end)=0; %q*sin(freqCarga*(t+dtdiv2)*2*pi);%20kN
    else
%         Rext(end)=0;
%         Rext2(end)=0;%20kN
    end
    %Integracion explicita
     Vnxt = V + dtdiv2*(M\(Rext(isFree)-C*V-K*D));
     Dnxt = D + dtdiv2*(V);
     % Integracion directa
%     k1V = M\(Rext(isFree)-C*V-K*D);
%     k1D = V;
%     
%     k2V = M\(Rext2(isFree)-C*(V+dtdiv2*k1V)-K*(D+dtdiv2*k1D));
%     k2D = V+dtdiv2*k1V;
%     
%     k3V = M\(Rext2(isFree) - C*(V+dtdiv2*k2V)-K*(D+dtdiv2*k2D));
%     k3D = V+dtdiv2*k2V;
%     
%     k4V = M\(Rext(isFree)-C*(V+dt*k3V)-K*(D+dt*k3D));
%     k4D = V+dt*k3V;
%     
%     Vnxt = V + dtdiv6*(k1V+2*(k2V+k3V)+k4V);
%     Dnxt = D + dtdiv6*(k1D+2*(k2D+k3D)+k4D);
    if sum(isnan(Dnxt))>0
        error('NaN found in solve')
    end 
    D_t(:,i+1) = D;
    D = Dnxt;
    V = Vnxt;
end

plot(D_t(end-1,:),'b.')

% pos = zeros(size(nodos));
% D = zeros(dof,1);
% for i=1:Nt
%     D(isFree) = D_t(:,i);
%     pos(:,1)=nodos(:,1)+D(1:2:end);
%     pos(:,2)=nodos(:,2)+D(2:2:end);
%     
%     scatter(pos(:,1),pos(:,2),'*k')
%     drawnow
%     
% end
