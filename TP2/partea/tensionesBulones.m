%% Analisis No lineal
% return
K=Kr;
M=Mr;
%% Non Linear analysis
Rext=zeros(dof,1);
Rext2=zeros(dof,1);
dofred=sum(isFree);
dt=.0001;
dtdiv2=dt/2;
dtdiv6=dt/6;
tfinal=1;
Nt=ceil(tfinal/dt);
applyDof=masa*6-3;
% dt = 1/max(omega)^2/4;
alpha = 0.002;
beta = 0.00002;%.1;

C = alpha*Mr+beta*Kr;
V = zeros(dofred,1); % Condicion inicial: velocidad = 0

D = zeros(dof,1); % desplazamientos = 0
%% Parametros de iteracion
iterskip=120;
D_t = zeros(dof,ceil(Nt/iterskip));

k=0;
for i = 0:Nt
    t = dt*i; %tiempo
    Rext(applyDof)=F0*sin(omegaexc*t);
    Rext2(applyDof)=F0*sin(omegaexc*(t+dt));
    %Integracion explicita
%      Vnxt = V + dtdiv2*(M\(Rext(isFree)-C*V-K*D(isFree)));
%      Dnxt = D(isFree) + dtdiv2*(V);
     % Integracion directa
    k1V = M\(Rext(isFree)-C*V-K*D);
    k1D = V;
    
    k2V = M\(Rext2(isFree)-C*(V+dtdiv2*k1V)-K*(D+dtdiv2*k1D));
    k2D = V+dtdiv2*k1V;
    
    k3V = M\(Rext2(isFree) - C*(V+dtdiv2*k2V)-K*(D+dtdiv2*k2D));
    k3D = V+dtdiv2*k2V;
    
    k4V = M\(Rext(isFree)-C*(V+dt*k3V)-K*(D+dt*k3D));
    k4D = V+dt*k3V;
    
    Vnxt = V + dtdiv6*(k1V+2*(k2V+k3V)+k4V);
    Dnxt = D + dtdiv6*(k1D+2*(k2D+k3D)+k4D);
    if sum(isnan(Dnxt))>0
        error('NaN found in solve')
    end 
    if mod(i,iterskip)==0
        k=k+1;
        D_t(:,k)=D;
    end
    D = Dnxt;
    V = Vnxt;
end
%% Grafico resultados
lims = [-.3, 1.3*max(nodos(:,1));
        -.3, 1.3*max(nodos(:,2));
        -2*max(nodos(:,3)) 2*max(nodos(:,3))];
pasaLaRepe=2;
mag=20;
Niter=size(D_t,2);
for j=1:pasaLaRepe
%     break %COMENTAR CUANDO SE QUIERE VER EL GRAFICO
for i=0:Niter-1
    t=i*iterskip*dt;
    posdef = nodos+mag*[D_t(1:6:end,i+1) D_t(2:6:end,i+1) D_t(3:6:end,i+1)];
    scatter3(posdef(:,1),posdef(:,2),posdef(:,3),'k.')
    hold on
    scatter3(lims(1,:),lims(2,:),lims(3,:),'w.')
    title('Analisis No lineal')
    ylabel(sprintf('Tiempo: %0.0f ms',t*1000))
    drawnow
    hold off
    pause(.01)
end
end

% Mis Vigas
E = 200e9;
rho = 7900;
nu=0.3;
A=b*h;
Iz = b*h^3/12;
Iy = h*b^3/12;
Jtors = h*b^3*(1/3-0.21*b/h*(1-b^4/12/h^4));
G=E/(2-2*nu);

SxBulones = zeros(size(rojos,1)*2+size(azules,1)*2,Niter);

% Parametros geometricos cuadrados (Cook pg 50. 2.9-6)
cy=1.5;cz=cy;cT=0.675*h;
clear x
xv = [1]; % puntos de interpolacion de tensiones. Solo un punto y es 
% al final de la viga porque ahi tengo mi bulon

basaMotor=[];
for r=reshape(rojos,[],1)' %BASAMENTO AL MOTOR
    aux=get_connectivity(r,elementos);
    if size(aux,1)>1
        basaMotor=[basaMotor;aux(1,:)];
    else
        basaMotor=[basaMotor;aux];
    end
end

fijoBasa=[];
for r=reshape(azules,[],1)' %Doble Fondo A basamento
    aux=get_connectivity(r,elementos);
    if size(aux,1)>1
        fijoBasa=[fijoBasa;aux(1,:)];
    else
        fijoBasa=[fijoBasa;aux];
    end
end

for i_t=1:Niter
    D=D_t(:,i_t);
    k=0;
for e=[basaMotor(:,1)',fijoBasa(:,1)']
    k=k+1;
    xv=1;
    n1 = nodos(elementos(e,1),:);
    n2 = nodos(elementos(e,2),:);
    [~, T]=vigorotar(zeros(12),n1,n2,[0 0 1]);
    Dlocal=T'*D(elemDof(e,:));
    % 1  2  3   4    5    6
    % u  v  w  phi  psi  tita
    u1=Dlocal(1);v1=Dlocal(2);w1=Dlocal(3);ph1=Dlocal(4); p1=Dlocal(5); t1=Dlocal(6);
    u2=Dlocal(7);v2=Dlocal(8);w2=Dlocal(9);ph2=Dlocal(10);p2=Dlocal(11);t2=Dlocal(12);


    Le=norm(n2-n1);
    ay=12*E*Iy/(Jtors*G*A*Le^2);
    az = 12*E*Iz/(Jtors*G*A*Le^2);
    by=1/(1-ay);
    bz = 1/(1-az);
    sig=0;
    for i=1:length(xv)
        x=xv(i);
        ddv = by*v1*(12*x - 6) - by*v2*(12*x - 6) - (by*t2*(ay - 6*x + 2))/2 + (by*t1*(ay + 6*x - 4))/2;
        ddw = bz*w1*(12*x - 6) - bz*w2*(12*x - 6) - (bz*p2*(az - 6*x + 2))/2 + (bz*p1*(az + 6*x - 4))/2;
%         dddv = 3*by*t1 + 3*by*t2 + 12*by*v1 - 12*by*v2;
%         dddw = 3*bz*p1 + 3*bz*p2 + 12*bz*w1 - 12*bz*w2;
        Nx = A*E/Le*(u2-u1); 
%         Tors = G*K*(ph2-ph1)/Le;
        My = E*Iy*ddw;
        Mz = E*Iz*ddv;
%         Vy = E*Iz*dddv;
%         Vz = E*Iy*dddw;
        % TENSIONES
        Sx = Nx/A-Mz*h/2/Iz-My*b/2/Iy;
%         Txy = Tors*cT/Jtors;
%         Ty = cy*Vy/A; Tz = cz*Vz/A;
        if k>length(basaMotor(:,1))
            SxBulones(k,i_t) = Sx;
        else
            SxBulones(k,i_t) = Sx;
        end
    end
end
end

maxSig=max(max(SxBulones)) %La tension maxima puede caer en un rigid link
maxSig2=max(max(SxBulones(length(basaMotor(:,1))+1:end,:))) %Tension NO rigid link
for k=1:size(SxBulones,1)
    for j=1:size(SxBulones,2)
        if SxBulones(k,j)==maxSig
            i=k;
            break
        end
    end
end
figure(2)
plot(dt*[0:Niter-1]*1000,SxBulones(i-1,:),'r^--')
title('Tensiones en el bulón más comprometido')
ylabel('\sigma_x [Pa]')
xlabel('Tiempo [ms]')

for k=1:size(SxBulones,1)
    for j=1:size(SxBulones,2)
        if SxBulones(k,j)==maxSig2
            i=k;
            break
        end
    end
end
figure(3)
plot(dt*[0:Niter-1]*1000,SxBulones(i-1,:),'r^--')
title('Tensiones en el bulón más comprometido no conectado a Rigid Beam')
ylabel('\sigma_x [Pa]')
xlabel('Tiempo [ms]')