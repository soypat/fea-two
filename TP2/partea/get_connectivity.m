function [connectedElem] = get_connectivity(nodo,elementos)
%
[Nelem, Nnodporelem]=size(elementos);
connectedElem=[];
for e=1:Nelem
    for n = 1:Nnodporelem
       if elementos(e,n)==nodo
           connectedElem=[connectedElem;e, n];
       end
    end
end

