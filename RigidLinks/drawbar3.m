function []=drawbar3(elementos,nodos)

Nelem=size(elementos,1);
Nnod = size(nodos,1);

scatter3(nodos(:,1),nodos(:,2),nodos(:,3),'r.');
hold on
for e = 1:Nelem
    n1=nodos(elementos(e,1),:);n2=nodos(elementos(e,2),:);
    plot3([n1(1) n2(1)],[n1(2) n2(2)],[n1(3) n2(3)],'k')
    hold on
end
hold off
end