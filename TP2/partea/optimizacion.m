%% Analisis No lineal y optimización
% return
K=Kr;
M=Mr;
%% Non Linear analysis
Rext=zeros(dof,1);
Rext2=zeros(dof,1);
dofred=sum(isFree);
dt=.0001;
dtdiv2=dt/2;
tfinal=6;
Nt=ceil(tfinal/dt);
applyDof=masa*6-3;
dt = 1/max(omega)/50000;
alpha = 0.02;
beta = 0.0002;%.1;

C = alpha*Mr+beta*Kr;
V = zeros(dofred,1); % Condicion inicial: velocidad = 0

D = zeros(dof,1); % desplazamientos = 0
%% Parametros de iteracion
iterskip=30;
D_t = zeros(dofred/2,ceil(Nt/iterskip));


for i = 0:Nt
    t = dt*i; %tiempo
    Rext(applyDof)=F0*sin(omegaexc*t);
    Rext2(applyDof)=F0*sin(omegaexc*(t+dt));
    %Integracion explicita
     Vnxt = V + dtdiv2*(M\(Rext(isFree)-C*V-K*D(isFree)));
     Dnxt = D(isFree) + dtdiv2*(V);
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
    if mod(i,iterskip)==0
        D_t(1:3:end,i+1)=D(1:6:end);
        D_t(2:3:end,i+1)=D(2:6:end);
        D_t(3:3:end,i+1)=D(3:6:end);
    end
    D = Dnxt;
    V = Vnxt;
end
% D_t(:,i+1) = 
iterskip=20;
% return
mag=50;
for i=0:iterskip:Nt
    t=i*iterskip*dt;
    posdef = nodos+mag*[D_t(1:3:end,i+1) D_t(2:3:end,i+1) D_t(3:3:end,i+1)]
    scatter3(posdef(:,1),posdef(:,2),posdef(:,3) )
    drawnow
    pause(.01)
end


