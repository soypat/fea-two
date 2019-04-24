function [w] = w_apoyado(x,y,a,b,N,P,Db)
% UGURAL 13.19 y 13.20
% Muy costoso numericamente, evitar uso
    w=0;
    syms xi eta
    pmn =@(m,n) eval(4/a/b * int( int( P*sin(m*pi*xi/a)*sin(n*pi*eta/b),xi,0,a ),eta,0,b));
    for m=1:2:N
        for n=1:2:N
           w=w+( pmn(m,n)/( (m/a)^2+(n/b)^2  )^2)*sin(m*pi*x/a)*sin(n*pi*y/b);
        end
    end
    w=w/(pi^4*Db);
end

