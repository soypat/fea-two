function [] = elPloto(nodos,elementos)
%ELPLOTO Plotea (nodos, elementos) con flechas orientadas
[Nelem, Nnodporelem]= size(elementos);  
% [Nnod, Ndim] = size(nodos);
xlim([min(nodos(:,1))-1 max(nodos(:,1))+1]);
ylim([min(nodos(:,2))-1 max(nodos(:,2))+1]);
scaleFact = 1.2;
positionFact = 1.6;
hold on
for e = 1:Nelem
    index=elementos(e,:);
    color = [randi(255) randi(255) randi(255)]/255;
    minX = min(nodos(:,1))*positionFact;
    minY = min(nodos(:,2))*positionFact;
    
    spanX = abs(max(nodos(:,1))-minX)*scaleFact;
    spanY = abs(max(nodos(:,2))-minY)*scaleFact;
   for n=1:Nnodporelem
      if n==1
          lastnod=index(Nnodporelem);
          
      end
      nod = index(n);
      eval(sprintf('ar%i%i = annotation("arrow");',e,n));
      hold on
      eval(sprintf('ar%i%i.X = [%f %f];',e,n,(nodos(lastnod,1)-minX)/spanX,(nodos(nod,1)-minX)/spanX));
      eval(sprintf('ar%i%i.Y = [%f %f];',e,n,(nodos(lastnod,2)-minY)/spanY,(nodos(nod,2)-minY)/spanY));
      eval(sprintf('ar%i%i.Color = [%f %f %f];',e,n,color(1),color(2),color(3)));
      hold on
      
%       arrow([nodos(lastnod,1) nodos(nod,1)],[nodos(lastnod,2) nodos(nod,2)],'Color',color);
      lastnod = nod;
      
   end
end
end

