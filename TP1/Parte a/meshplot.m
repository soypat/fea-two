function meshplot(elementos,nodos,color,verbose)
% MESHPLOT  Graficador de mallas.
%
% MESHPLOT(elementos,nodos,color)
%
% Numeración de los nodos de los elementos:
%  4---7---3
%  |       |
%  8   9   6
%  |       |
%  1---5---2
%
% elementos: matriz de conectividades.
% nodos:     matriz de coordenadas nodales.
% color:     string que especifica el color para los bordes de los
%            elementos.

nNodos = size(elementos,2);
switch nNodos
    case {8,9}
        vNod = [1 5 2 6 3 7 4 8];
    otherwise
        vNod = 1:4;
end

h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1])
daspect([1 1 1])

%% Display numeracion nodos y elementos
if verbose
    for n=1:size(nodos,1)
        text(nodos(n,1),nodos(n,2),num2str(n))
    end
    
    for e=1:size(elementos,1)
        text(mean(nodos(elementos(e,:),1))-5,mean(nodos(elementos(e,:),2)),num2str(e),'Color','r')
    end
end
