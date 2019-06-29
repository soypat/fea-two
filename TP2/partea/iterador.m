%% PROBLEM
LargoMotor = 5.719;
% Modifiables
AB = [1.8 1.8];
longE = LargoMotor/4*ones(3,1)';%[1.2 2.95 ] %E1, E2, E3
rpmexc = 600;
omegaexc = rpmexc/60*2*pi;
Lelemax=.125; %Longitud máxima de elementos

hv = 0.070:.00125:.1;
bv = 0.045:.00125/2:0.06;
Nh = length(hv);
Nb = length(bv);

Niter=Nh*Nb;
ih = 1;
ib = 1;
AModhb = zeros(Nh,Nb);
Shb = zeros(Nh,Nb);
AProphb = zeros(Nh,Nb);
Secchb = zeros(Nh,Nb);

l = waitbar(0,'Beginning iteration');
Ntotal = Nh*Nb;

for h = hv 
    for b = bv
        waitbar(((ih-1)*Nb+ib)/Ntotal,l,sprintf('h=%0.3f[m]',h))
        mainiterable
%         Shb(ih,ib) = maxsig;
        AModhb(ih,ib) = AmpMod;
        AProphb(ih,ib) = AmpProp;
        Secchb(ih,ib) = h*b;
        ib=ib+1;
    end
    ib=1;
    ih=ih+1;
end
% close(l)
close all
[Xb ,Yh] = meshgrid(bv,hv);
figure(1)
surf(Xb,Yh,AModhb)
title('Amplitud en funcion de h y b (Amort. modal)')
xlabel('b [m]')
ylabel('h [m]')

figure(2)
isnn=AProphb==0;
AProphb(isnn)=nan;
surf(Xb,Yh,AProphb)
title('Amplitud en funcion de h y b (Amort. prop.)')
xlabel('b [m]')
ylabel('h [m]')

figure(3)
surf(Xb,Yh,AModhb./Secchb)
title('Seccion en funcion de h y b')
xlabel('b [m]')
ylabel('h [m]')

%% Find minimum amplitud
mini=min(min(AProphb));
index=find(AProphb==mini);
ib = ceil(index/Nb);
index=find(AProphb'==mini);
ih = ceil(index/Nh);

hopt = hv(ih);
bopt = bv(ib);

fprintf('h optimo=%0.0fmm\nb optimo=%0.0fmm',hopt*1000,bopt*1000) % No significa que sea la mejor opción
