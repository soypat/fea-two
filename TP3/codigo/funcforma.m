%% Func Forma Hexahedral 8 node elements

syms ksi eta zeta x y z
unod =[...
    -1    -1    -1
    -1     1    -1
     1     1    -1
     1    -1    -1
    -1    -1     1
    -1     1     1
     1     1     1
     1    -1     1];
 X = [1 x y z x*y x*z y*z x*y*z];
 
A=nan(length(X));
k=0;
for inod=unod'
   k=k+1;
   x=inod(1);y=inod(2);z=inod(3);
   A(k,:)=subs(X);
end
x=ksi;
y=eta;
z=zeta;
X=subs(X);
N = X/A;
NL(1,1:3:3*length(N))=N;
NL(2,2:3:3*length(N))=N; %Tiene la forma de las funciones de forma encontradas en el cook pg 206, ecuacion (6.2-2). 
NL(3,3:3:3*length(N))=N;
dN = [diff(NL(1,:),ksi); diff(NL(2,:),eta); diff(NL(3,:),zeta)];
dNaux=[diff(N,ksi);diff(N,eta);diff(N,zeta)]; %Para calcular jacobiano
