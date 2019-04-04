%% SHape Funspara Placa Kirchoff C1

clear dN N dNaux X dNxx dNyy dNxy 
syms x y real

X = [1 x y x.^2 x.*y y.^2 x.^3 (x.^2).*y x.*(y.^2) y.^3 x^3*y x.*(y^3)];
% w1 t1 t2 
Xdx = diff(X,x);
Xdy = diff(X,y);
uDof = [-1 -1
        -1 -1
        -1 -1
        1 -1
        1 -1
        1  -1
        1  1
        1  1
        1  1
        -1 1
        -1 1
        -1 1];
uNod =[-1 -1;1 -1;1 1;-1 1];
A = zeros(size(uDof,1),length(X));

for i=1:size(uDof,1)
    x=uDof(i,1); y = uDof(i,2);
    if mod(i,3)==0
        A(i,:) = double(subs(Xdy));
    elseif mod(i,3)==2
        A(i,:) = double(subs(Xdx));
    elseif mod(i,3)==1
        A(i,:) = double(subs(X));
    else
        error('WANRASDKDA')
    end
    
end

syms x y real
shapefuns = X*inv(A);
dNxx = diff(shapefuns,x,x);
dNyy = diff(shapefuns,y,y);
dNxy = diff(shapefuns,x,y);

B = [dNxx; dNyy; 2*dNxy];
N = shapefuns;
% for n = uNod'
%    x=n(1);y=n(2);
%    n
%    sum(subs(N));
%    subs(N)
% end
