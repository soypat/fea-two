function meshplot(elementos,nodos,bc,color,verbatim)
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

L = abs(max(nodos(:,1))-min(nodos(:,1)));

nNodos = size(elementos,2);
switch nNodos
    case {8,9}
        vNod = [1 5 2 6 3 7 4 8];
    otherwise
        vNod = 1:nNodos;
end

h1 = patch('Faces',elementos(:,vNod),'Vertices',nodos);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'XColor',[1 1 1],'YColor',[1 1 1])
daspect([1 1 1])
hold on
for n=1:size(nodos,1)
    if bc(n,1)==1
        fill([nodos(n,1);nodos(n,1)-0.06*L;nodos(n,1)-0.06*L],[nodos(n,2); nodos(n,2)-0.03*L; nodos(n,2)+0.03*L],[1;1;1])
    end
    if bc(n,2)==1
        fill([nodos(n,1);nodos(n,1)-0.03*L;nodos(n,1)+0.03*L],[nodos(n,2); nodos(n,2)-0.06*L; nodos(n,2)-0.06*L],[1;1;1])
    end
    if verbatim
        text(nodos(n,1),nodos(n,2),num2str(n),'FontSize',12,'Color','r')
    end
end
if verbatim
    for e=1:size(elementos,1)
        text(mean(nodos(elementos(e,:),1)),mean(nodos(elementos(e,:),2)),num2str(e),'FontSize',12,'Color','b')
    end
end