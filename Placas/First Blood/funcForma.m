%% SHape Funspara Placa Kirchoff C1

clear dN N dNaux X dNxx dNyy dNxy 
syms x y real

X = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3];
% w1 t1 t2 
Xdx = diff(X,x);
Xdy = diff(X,y);
uNod = [-1 -1
        -1 -1
        -1 -1
        -1 1
        -1 1
       -1  1
        1  1
        1  1
        1  1
        1 -1
        1 -1
        1 -1];
A = zeros(size(uNod,1),length(X));

for i=1:size(uNod,1)
    x=uNod(i,1); y = uNod(i,2);
    if mod(i,3)==0
        A(i,:) = double(subs(Xdy));
    elseif mod(i,2)==0
        A(i,:) = double(subs(Xdx));
    else
        A(i,:) = double(subs(X));
    end
   
end

syms x y real
shapefuns = X/A;
dNxx = diff(shapefuns,x,x);
dNyy = diff(shapefuns,y,y);
dNxy = diff(shapefuns,x,y);

B =[dNxx; dNyy; 2*dNxy];
