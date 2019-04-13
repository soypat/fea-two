function [w] = w_analytic(x,y,a,b,N,P,Db)
    w=0;
    for m=1:2:N
        for n=1:2:N
            w=w+16*P/(pi()^6*Db).*(sin(m.*pi().*x./a).*sin(n.*pi().*y./b))./(m*n*((m/a)^2+(n/b)^2)^2);    
        end
    end
end