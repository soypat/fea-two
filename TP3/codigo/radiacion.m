
%% Iteracionar al radiacion
Rrad=zeros(dof,1);
for e = 1:Nelem
    index = elementos(e,:);
    elenod = nodos(index,:);
    s=0;
    for snod = supnod'
        s=s+1;
        xnod=elenod(snod,1); ynod=elenod(snod,2); znod=elenod(snod,3);
        if sum(  xnod==0 | xnod==0.8 | ynod==0 | znod ==.8 | znod==0  )>2
            k=k+1;
            r = zeros(4,1);
            supindex = index(supnod(s,:));
            Tupg=0;
            for ipg = 1:npgs %Acá comienza la integracion sobre superficie
                J = dNauxs{s}{ipg}*elenod;
                % con la optimizacion se desvirtua todo, pero bue, es lo
                % que hay
                for n=1:Nnodporelem
                    Tupg = Tupg + T(index(n))*Ns{s}(ipg,n);
                end
%                 Tupg = T(index)'*Ns{s}(ipg,:)'; %interpolacion
                r = r-Ns{s}(ipg,supnod(s,:))'*boltz*(Tupg^4 - Trad)*det(J)*wpg(ipg);%Ojo, es negativo porque el calor se ``va''
            end
            Rrad(supindex)=Rrad(supindex)+r;
        end
    end
end