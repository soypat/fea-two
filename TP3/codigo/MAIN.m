clear N dN dNaux NL ksi eta zeta x y z wn wz wy wx xv xnv
% clc
%%PROBLEM
Q = 2e3; %2kW/m^3 generados
L=0.8;
volumen = L*L*L/8;

%% SOLVE:
funcforma

% MODIFICAS CANTIDAD DE ELEMENTOS
div = 19;
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
     0 0 k]; %considero material isotropico 3D

%% Punto Gauss con método acelerado de calculo
[wpg, upg, npg] = gauss([2 2 2]); % Ver tabla 6.8-1 Cook
[Ns,~,dNauxs,~] = shapefunGP(upg,N,dN,dNaux,NL); %Funcion optimizadora

%% Rigidity Assemble and charges
Kg = zeros(dof,dof);
Cg = zeros(dof,dof);
Rgen=zeros(dof,1);
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
    Rgen(meindof) = Rgen(meindof)+r; %cargas térmicas por generacion
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
K = Kg;
C= Cg;
R = zeros(dof,1);
Rk = zeros(dof,1);
T = zeros(dof,1);

%% Condiciones de borde
cc=false(dof,1);
cc(medio)=true;
xx = ~cc;
T(medio) = 300; %K condicion de borde
Tlast = T;

%% LAST RESOLUTION
% v = load('Tf.mat');
% T = v.T;
%% Iteracion para llegar a regimen estacionario (estado estable)
prepRad %Preparo optimizador de radiación (guardo funcformas en puntos gauss)
keepGoing=true;
i=1;
T(xx)=2.7;%Condicion inicial para que converga
relajacion = 8;
relajar=1/relajacion;

while keepGoing
    %% SOLVER RADIACIÓN
%     convectar %Se genera el vector {Rrad} en funcion de {T}
    radiacion
    R = Rgen + Rrad; 
    T(xx) = K(xx,xx)\(R(xx)-K(xx,cc)*T(cc));
    DT = relajar*(T-Tlast);
    T=Tlast+DT;
    
%     AreaRad % Area radiativa es 0.48m^2 ... verificar variable
    %% PostProcess (graficado)
    if norm(Tlast-T)/norm(Tlast) <1e-7 && i>20
        warning('Ending simulation. error low')
        keepGoing=false;
    elseif mod(i,1)==0
        Tin = T(interiornod);
        scatter(i,Tin(end),'r.')
        hold on
        scatter(i,T(1),'b.')
        hold on
        drawnow
    end
    Tlast = T;
    if isnan(T(wallnod(1)))
        error('NAN FOUND!')
        break
    end
    i=i+1;
end

title(sprintf('Convergencia con %0.0f elementos',Nelem))
ylabel('Temperatura [K]')
xlabel('Iteracion')
legend('T_i','T_s')
grid on
figure(2)
plot(xv,T(xnv),'k--^')
title('Perfil de temperaturas Centro-Superficie')
xlabel('x [m]')
ylabel('Temperatura [K]')
grid on

figure(3)
plot(wz,T(wn),'m--s')
title('Perfil de temperaturas a lo largo de superficie')
ylabel('Temperatura [K]')
xlabel('Posicion sobre superficie [m]')
grid on
figure(4)
plot(wz,boltz*(T(wn)-2.7^4),'b--x')
title('Calor intercambiado de la superficie con el espacio')
ylabel('Calor radiado por unidad de area [W/m^2]')
xlabel('Posicion sobre superficie [m]')
grid on

save('Tf.mat','T')

%% Aplicamos condiciones de borde conocidas 
%  para obtener un pantallazo del problema
return
title(sprintf('Oscilación con %0.0f elementos sin calor generado',Nelem))
ylabel('Temperatura [K]')
xlabel('Iteracion')
legend('T_i','T_s')
grid on


title(sprintf('Divergencia con %0.0f elementos',Nelem))
ylabel('Temperatura [K]')
xlabel('Iteracion')
legend('T_i','T_s')
grid on
%% Reso Radiacion


%% Graph problema

% scatter3(nodos(:,1),nodos(:,2),nodos(:,3),'bo')
% hold on
% scatter3(nodos(medio,1),nodos(medio,2),nodos(medio,3),'r*')
% hold on
% scatter3(nodos(wallnod,1),nodos(wallnod,2),nodos(wallnod,3),'m.')
% hold on
% scatter3(nodos(interiornod,1),nodos(interiornod,2),nodos(interiornod,3),'g^')