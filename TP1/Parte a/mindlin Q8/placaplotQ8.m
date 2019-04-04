function [] = placaplotQ8(nodos,D)
% Ploteate los Q8 como un campeon
Nnod = size(nodos,1);
w = D(1:3:end);
if length(w)~= Nnod
   error('Inconsistencia dimensiones nodos y desplazamientos en placaplotQ8')
end
scatter3(nodos(:,1)',nodos(:,2)',w)

end

