A = Mr\Kr; 

[Vr, eigVal]=eig(A);

V=zeros(dof);

V(isFree,isFree)=Vr;

Db = Vr;

dofred = size(Vr,2); % Lo mismo que %dofred=sum(isFree)

Phi = zeros(dofred,dofred);
% omega = zeros(dofred,1);
omega = sqrt(diag(eigVal));
for i = 1:dofred
    Dbi=Db(:,i);
    
    aux=Dbi'*Mr*Dbi;
    
    Phi(:,i) = Dbi/sqrt(aux);
    Phii = Phi(:,i); %Phii = []
%     omega(i) = sqrt((Dbi' * Kr *Dbi)/aux);% Cook (11.4-13)
end
% ESP  = Phi' * Kr *Phi;
% omega2 = sqrt(diag(ESP)); % y una cuarta para que tengas

Nmodos = 3;
omega2 = omegaexc*1.15;
omega1 = omega(end);
ksi1=0.073;
ksi2=0.2;
%% Amortiguamiento Proporcional
alpha = 2*omega1*omega2*(ksi1*omega2-ksi2*omega1)/(omega2^2 - omega1^2);
beta = 2*(ksi2*omega2-ksi1*omega1)/(omega2^2 - omega1^2);
Cprop = diag(alpha*eye(dofred,dofred)+beta*eigVal);
%% Amortiguamiento Modal
ksiModal = 0.02;
dampingModal = ksiModal*ones(dofred,1);
Cmodal = (2*dampingModal'*eigVal)';

%% Obtencion de amplitudes de oscilacion
Rmod = -Phi'* Rrexc;
chi = omegaexc./omega;

Zmod = (Rmod./omega.^2)./ ...
                sqrt((1-chi.^2).^2+(2.*Cmodal.*chi).^2);
Zprop = (Rmod./omega.^2)./ ...
                sqrt((1-chi.^2).^2+(2.*Cprop.*chi).^2);

AmpMod = sum(abs(Zmod));
AmpProp = sum(abs(Zprop));

return
% input_omega = 0:2000:200000;
input_omega = 0:1:(omegaexc*1.25);
Nfrec = length(input_omega);
Amod = zeros(Nfrec,1);
Aprop = zeros(Nfrec,1); 


for f = 1:Nfrec   
        chi = input_omega(f)/omega(i);
        Zmod = (Rmod./omega.^2)./ ...
            sqrt((1-chi.^2).^2+(2.*Cmodal.*chi).^2);
        Zprop = (Rmod./omega.^2)./ ...
            sqrt((1-chi.^2).^2+(2.*Cprop.*chi).^2);

        Amod(f) = abs(sum(Zmod));% tomo valor absoluto para que me quede lindo
        Aprop(f) = abs(sum(Zprop));
end
figure(3)
oe = omegaexc/2/pi*[1 1];
am = [min(min(Amod,Aprop));max(max(Amod,Aprop))];
semilogy(oe,am)
hold on
semilogy(input_omega/(2*pi),Amod)
hold on
semilogy(input_omega/(2*pi),Aprop)

grid on
xlabel('Frecuencia excitacion [Hz]')
ylabel('Amplitud [UA]')

hold on

legend('\omega_e','Modal','Proporcional')
