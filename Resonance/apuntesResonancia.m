% Apuntes

%% Reducción de orden n x n ---> n x m ( n del cook es mi dofred)
% m_red=6; %Reduzco mi problema a m eigenvalores
% for i = dofred-1:-1:dofred-1-m_red
%     Di = D
% end

%% Obtencion de desplazamientos modales

% Z = Phi' *Mg(isFree,isFree) ;



%% CARGAS OTRA FORMA:
% El primer modo se lleva Rmodal(1) cantidad de carga 
% Ver descomposicion modal - cargas armonicas
% Lo que estamos haciendo es armar la maquina que va hacer un sine sweep
% al modelo que armamos,

% Si lo hago de la otra forma 1) v= vbar sin wt:

% [Mcc  Mcx]*{D''c} + [Ccc  Ccx]*{D'c} + [Kcc ...]*Dc = 0
% [Mxc  Mxx] {D''x}   [Cxc  Cxx] {D'x}             Dx      % Descomposicion
% conocido/desconocido es una idea loca. hay forma mas facil:
% MxxD''x+CxxD'x + KxxDx = -MxcD''c -CxcD'c -KxcDc
%                        =+Mxcw^2Asen(wt) - CxcAwcos(wt)-KxcAsen(wt)
