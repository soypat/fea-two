function [nodos,elementos] = mesh1D(nodos,n1,n2,Nele)
%MESHV mesh beween nodes
Nnodporelem = 2;
elementos=zeros(Nele,Nnodporelem);
if Nele == 1
   elementos(1,:) = [n1 n2];
   return
else
    Nnod = size(nodos,1);
    N1 = nodos(n1,:);
    N2 = nodos(n2,:);
    div = Nele-1;
    newnodos = zeros(div,3);  
    xv = linspace(N1(1),N2(1),Nele+1);
    yv = linspace(N1(2),N2(2),Nele+1);
    zv = linspace(N1(3),N2(3),Nele+1);
    
    for i=1:Nele
        
        if i==1
            newnodos(i,:) = [xv(i+1) yv(i+1) zv(i+1)];
            elementos(i,:) = [n1 Nnod+i];
        elseif i==Nele
            elementos(i,:) = [Nnod+i-1 n2];
            break
        else
            newnodos(i,:) = [xv(i+1) yv(i+1) zv(i+1)];
            elementos(i,:) = [Nnod+i-1 Nnod+i];
        end
    end
    nodos = [nodos;newnodos];
end
end

