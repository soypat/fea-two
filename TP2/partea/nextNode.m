function [type] = nextNode(x,X1,X2,X3)
k1=1;
k2=1;
k3=1;
N1= length(X1);
N2= length(X2);
N3= length(X3);
while X1(k1)<=x
    k1=k1+1;
    if k1>N1
        k1=N1; X1(N1)=inf;
    end
end

while X2(k2)<=x
    k2=k2+1;
    if k2>N2
        k2=N2; X2(N2)=inf;
    end
end

while X3(k3)<=x
    k3=k3+1;
    if k3>N3
        k3=N3; X3(N3)=inf;
    end
end

if X1(k1)<X2(k2) && X1(k1)<X3(k3)
    type=1;
    return
elseif X1(k1)>X2(k2) && X2(k2)<X3(k3)
    type=2;
    return
elseif X1(k1)>X3(k3) && X2(k2)>X3(k3)
    type=3;
    return
else 
    error('coincident nodes')
end
end