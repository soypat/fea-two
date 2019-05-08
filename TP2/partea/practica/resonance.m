A = Mr\Kr; 

[Vr, eigVal]=eig(A);

V=zeros(dof);

V(isFree,isFree)=Vr;

Db = Vr;

dofred = size(Vr,2); % Lo mismo que %dofred=sum(isFree)

Phi = zeros(dofred,dofred);
omega = zeros(dofred,1);

for i = 1:dofred
    Dbi=Db(:,i);
    
    aux=Dbi'*Mr*Dbi;
    
    Phi(:,i) = Dbi/sqrt(aux);
    Phii = Phi(:,i); %Phii = []
    omega(i) = sqrt((Dbi' * Kg(isFree,isFree) *Dbi)/aux);% Cook (11.4-13)
end
ESP  = Phi' * Kg(isFree,isFree) *Phi;
omega2 = sqrt(diag(ESP)); % y una cuarta para que tengas
