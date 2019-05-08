
Ndofpornod = 2;
Nelem = 10;
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
A = Mg(isFree,isFree)\Kg(isFree,isFree); %F de frecuencia


Kr = Kg(isFree,isFree);
Mr = Mg(isFree,isFree);

[Vr, eigVal]=eig(A);

V=zeros(dof);
V(isFree,isFree)=Vr;


Db = Vr;
%Cantidad de dof despues de reduccion
dofred = size(Vr,2); % Lo mismo que %dofred=sum(isFree)

% PLOTEO ALGO
Ug = flip(flip(V,1),2); % Desplazamientos globales Ug (incluye angulos)
Uy = Ug(1:2:end,:); %Desplazamientos en y 

n = 4; % Numero de modo.
plot(nodos(:,1),Uy(:,n));



Phi = zeros(dofred,dofred);
omega = zeros(dofred,1);
omegaray = zeros(dofred,1);
ray = zeros(dofred,1);
for i = 1:dofred
    Dbi=Db(:,i);
    
    aux=Dbi'*Mr*Dbi;
    
    Phi(:,i) = Dbi/sqrt(aux);
    Phii = Phi(:,i); %Phii = []
    
    omegaray(i) = sqrt((Dbi' * Kg(isFree,isFree) *Dbi)/aux);% Cook (11.4-13)
    ray(i)= sqrt((Phii' * Kg(isFree,isFree) *Phii)/(Phii' * Mg(isFree,isFree) *Phii)); %Cook (11.7-1)b
    omega(i) = sqrt(Phii'* Kg(isFree,isFree)*Phii); % Idem
    % Wow son identicas las tres! como es? magia?  <No boludo, es matematica. si M=1 entonces obvio>
end
ESP  = Phi' * Kg(isFree,isFree) *Phi;
omega2 = sqrt(diag(ESP)); % y una cuarta para que tengas

%[omegaray ray omega omega2] %Descomentar para ver que son identicas

%% Cargas Externas / cargas modales (P y Rmodal [son lo mismo])
input_omega = 1:1:10000;
Nfrec = length(input_omega);

acel = @(x,A,omega) A.*sin(omega.*x);
Rext = zeros(dof,1); % (11.7-5) cook
Rmodalenfrecuencia = zeros(dofred,Nfrec); % (11.7-5) cook
for f =1:Nfrec
    for i=1:2:dof % los dof en y los sacudo
        n=round((i+1)/2);
        x = nodos(n,1);
        Rext(i) = 1;%;acel(x,1,input_omega(f));
%         if abs(x-L)<1e-4
%             Rext(i) = 1;%acel(x,1,input_omega(f));
%         end
    end
    Rextred = Rext(isFree);
    Rmodalenfrecuencia(:,f) = -Phi'* Mg(isFree,isFree)*Rextred;
end
% Rext= - Mg(isFree,isFree)*Rext;
% Rmodal = Phi' *Rext;  %Cada modo se lleva una carga. asi funciona esta wea
% OTRA FORMA:

%Comienza la magia negra
  % sine sweep dominio frecuencias
input_ksi = 0.02:.05:.2; %damping
Nksi = length(input_ksi);

Amp = zeros(Nfrec,Nksi);

for k = 1:Nksi
    for f = 1:Nfrec
        z=zeros(dofred,1);
        for i=1:dofred
            x = floor( (i+2)/2 );
            beta = input_omega(f)/omega(i);
            z(i) = (Rmodalenfrecuencia(i,f)/omega(i)^2)/sqrt((1-beta^2)^2+(2*input_ksi(k)*beta)^2);
        end
        Amp(f,k) = abs(sum(z));% tomo valor absoluto para que me quede lindo
    end
    semilogy(input_omega,Amp(:,k))
    hold on
end

%% Damping modal (Ni lo uso)
% La matriz de damping te queda diagonal... Genial! Atentos a la ecuacion
% (11.5-3) del cook
xi1= 0.02; %Damping lo escojo como 0.02 para ambas frecuencias que se yo.
xi2= xi1;
w1 = omega(end);
w2 = omega(end-1);
alfa = 2*w1*w2*(xi1*w2-xi2*w1)/(w2^2-w1^2);
beta = 2*(xi2*w2-xi1*w1)/(w2^2-w1^2);

Cmodal = Phi'*(alfa*Mg(isFree,isFree) +beta*Kg(isFree,isFree))*Phi;





