function [w] = w_max(a,b,N,P,Db)
%   SOLO PARA CARGA UNIFORME APOYADA
    w=0;
    for m=1:2:N
        for n=1:2:N
            w=w+((-1)^(m/2+n/2-1)/( (m/a)^2+(n/b)^2  )^2)/m/n;
        end
    end
    w=16*P*w/(pi^6*Db);
end
