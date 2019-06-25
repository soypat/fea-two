clear N dN dNaux NL ksi eta zeta x y z wn wz wy wx xv xnv
clc
%%PROBLEM
Q = 2e3; %2kW/m^3 generados
L=0.8;
volumen = L*L*L/8;

%% SOLVE:
funcforma

div = 11;
div2=ceil(9/2);

[nodos, ~,elementos] = mesh3D([0 L/2;0 L/2;0 L/2],[div div div]);
medio = size(nodos,1);%ceil(div^3/2);


Ndofpornod = 1; % T
[Nnod, Ndim]     =     size(nodos);
[Nelem, Nnodporelem] = size(elementos);
Ndofporelem  =  Ndofpornod*Nnodporelem;
dof = Ndofpornod*Nnod;
DOF = reshape(1:dof,Ndofpornod,[])';

%% Material
cp=528; %J/kg.K
rho =4500;
cap = cp*rho;
capTot = cap*volumen; % deberia dar igual que la suma de toda mi matriz capacidad
% fijate sum(sum(C))
k=17;%W/m/K
Kc = [k 0 0 
     0 k 0
     0 0 k]; %considero material isotropico

%% Punto Gauss con método acelerado de calculo
[wpg, upg, npg] = gauss([2 2 2]); % Ver tabla 6.8-1 Cook
[Ns,~,dNauxs,~] = shapefunGP(upg,N,dN,dNaux,NL);

%% Rigidity Assemble and charges
Kg = zeros(dof,dof);
Cg = zeros(dof,dof);
R=zeros(dof,1);
wallnod=false(Nnod,1);
for e = 1:Nelem
    Ke = zeros(Ndofporelem);
    Ce = zeros(Ndofporelem);
    r=zeros(Ndofporelem,1);
    index = elementos(e,:);
    elenod = nodos(index,:);
    for ipg = 1:npg
        J    = dNauxs{ipg}*elenod;
        dNxyz = J\dNauxs{ipg};
        B = [dNxyz(1,:);dNxyz(2,:);dNxyz(3,:)];
        
        Ce = Ce + Ns{ipg}'*cap*Ns{ipg}*det(J)*wpg(ipg);
        Ke = Ke + B'*Kc*B*wpg(ipg)*det(J);
        r = r+Ns{ipg}'*Q*wpg(ipg)*det(J);
    end
    meindof = reshape(DOF(index,:)',1,[]);
    R(meindof) = R(meindof)+r; %cargas térmicas por generacion
    Cg(meindof,meindof)=Cg(meindof,meindof)+Ce;
    Kg(meindof,meindof)=Kg(meindof,meindof)+Ke;
end

k=1;
wk=1;
interiornod = false(Nnod,1);
for n=1:Nnod
    x=nodos(n,1);y=nodos(n,2);z=nodos(n,3);
    if x==0 || z==0 || y==0
        wallnod(n)=true;
    elseif n~=Nnod
        interiornod(n)=true;
    end
    if x>=0 && y==L/2 && z==L/2
        xv(k) = x;
        xnv(k)=n;
        k=k+1;
    end
    if x==0 && y==L/2
        wn(wk)=n;
        wy(wk)=y;
        wz(wk)=z;
        wx(wk)=x;
        wk=wk+1;
    end
end
intnod= find(interiornod);

Rgen = R;
K = Kg;
C= Cg;
%% Resuelvo Condiciones iniciales
Ts = 2.7;
T = zeros(dof,1);
cc = wallnod;
% T(interiornod)=(Ts+300)/2;
T(cc) = Ts;
cc(medio)=true;
T(medio) = 300; %K condicion de borde
xx = ~cc;
Tx = K(xx,xx)\(R(xx)-K(xx,cc)*T(cc));
Qx = K(cc,xx)*Tx+K(cc,cc)*T(cc);
R(cc)=Qx;
T(xx)=Tx;

%% LAST RESOLUTION
% v = load('Tf.mat');
% T = v.T;
%% Transitorio
dt = 150;
t_tot = 100000;
Nt = ceil(t_tot/dt);
prepRad
dtdiv2=dt/2;
beta=0.5;
keepGoing=true;
i=1;
iterskip=30;
while keepGoing
    tiempo = i*dt;
    radiacion
    Rnxt = Rgen+Rrad;
    if beta == 0
        Tnxt = C\((C-dt*K)*T + dt*R );
    elseif beta==.5
        Tnxt = (C+dtdiv2*K)\((C-dtdiv2*K)*T + dtdiv2*(R+Rnxt ));
    else
        Tnxt = (C+beta*dt*K)\((C-dt*(1-beta)*K)*T + dt*((1-beta)*R+beta*Rnxt));
    end
%     norm(T-Tnxt)/norm(T)
    if norm(T-Tnxt)/norm(T) <1e-8 && i>20
        warning('Ending simulation. error low')
        keepGoing=false;
    elseif mod(i,iterskip)==0
        scatter(dt*i,T(intnod(end)),'r.')
        hold on
        scatter(dt*i,T(1),'b.')
        hold on
        drawnow
    end
    
    T = Tnxt;
    T(medio)=300;
    if isnan(T(62))
        error('NAN FOUND!')
        break
    end
    R = Rnxt;

%     drawnow
    i=i+1;
end
return
title(sprintf('Convergencia con %0.0f elementos',Nelem))
ylabel('Temperatura [K]')
xlabel('Tiempo [s]')
legend('T_i','T_s')
grid on

plot(xv,T(xnv),'k--^')
title('Perfil de temperaturas Centro-Superficie')
xlabel('x [m]')
ylabel('Temperatura [K]')
grid on


plot(wz,T(wn),'m--s')
title('Perfil de temperaturas a lo largo de superficie')
ylabel('Temperatura [K]')
xlabel('Posicion sobre superficie [m]')
grid on
plot(wz,boltz*(T(wn)-2.7^4),'b--x')
title('Calor intercambiado de la superficie con el espacio')
ylabel('Calor radiado por unidad de area [W/m^2]')
xlabel('Posicion sobre superficie [m]')
grid on

save('Tf.mat','T')

%% Aplicamos condiciones de borde conocidas 
%  para obtener un pantallazo del problema





%% Reso Radiacion


%% Graph problema

% scatter3(nodos(:,1),nodos(:,2),nodos(:,3),'bo')
% hold on
% scatter3(nodos(medio,1),nodos(medio,2),nodos(medio,3),'r*')
% hold on
% scatter3(nodos(wallnod,1),nodos(wallnod,2),nodos(wallnod,3),'m.')
% hold on
% scatter3(nodos(interiornod,1),nodos(interiornod,2),nodos(interiornod,3),'g^')