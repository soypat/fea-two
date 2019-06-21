%% CARGAS
R=zeros(dof,1);


supnodos=[1 2 3 4; %superficie 1
          1 2 6 5; %Estos son los nodos que conforman geometria ladrillo/hexahedro
          2 3 7 6;
          3 4 8 7;
          1 4 8 5;
          5 6 7 8]; %superficie 6
 supnodosextended=[1 2 3 4 9 10 11 12; %superficie 1
          1 2 6 5 9 18 13 17;
          2 3 7 6 10 19 14 18;
          3 4 8 7 11 20 15 19;
          1 4 8 5 12 20 16 17;
          5 6 7 8 13 14 15 16]; %superficie 6
npg=9;% GAUSS PARA SUPERFICIES n=2
wpg=reshape([5/9;8/9;5/9]*[5/9 8/9 5/9],1,[]);

sup=0;
olds=0;
r_time=tic;
Area=0;
k=0;
for e = 1:Nelem
    index = elementos(e,:);
    elenod = nodos(index,:);
    
    s=0;
    
    for snod = supnodos'
        s=s+1;
        xnod=elenod(snod,1);      % USAR ynod, znod etc.
        if sum(abs(   xnod-16   )<1e-4)>2 %MODIFICAR ESTA CONDICION PARA SUPERFICIE
            if olds~=s %Esto es un optimizador, no hace falta guardar que superficie ando
                olds=s;
                [surf_upg, VV] = getsurfupg(s);
                [Ns, ~ ,dNauxs, ~ ] = shapefunGP(surf_upg,N,dN,dNaux,NL);
                break
            else
                break
            end
        else
            if s==6
                s=7; %Llegamos a la ultima superf. y no encontramos aplicacion de carga
                continue
            end
        end
    end
    if s==7;continue;end %Si no encontro superficie, chau elemento.
    meindof = reshape(DOF(index,:)',1,[]);
    k=k+1;
    for ipg = 1:npg
        J    = dNauxs{ipg}*elenod;
        Area=Area+det(J)*wpg(ipg);
        VxV = cross(J(VV(1),:),J(VV(2),:)).';%Vector Columna. Ver Cook 6.9-8
        for snod = supnodosextended(s,:)
            Rindex = [snod*3-2 snod*3-1 snod*3];
            R(meindof(Rindex)) = R(meindof(Rindex)) ...
                +Ns{ipg}(snod)*qx*VxV*wpg(ipg);
        end
    end
end
fprintf('\nCargas obtenidas en %0.2f segundos\n\n',toc(r_time))
defuzz_time=tic;
R=defuzz(R);
toc(defuzz_time)