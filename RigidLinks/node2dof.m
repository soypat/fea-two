function elemDof = node2dof(elementos,Ndofpornod)
  [Nelem, Nnodporelem] = size(elementos);
  elemDof = zeros(Nelem,Nnodporelem*Ndofpornod);
  for e=1:Nelem
    for n=1:Nnodporelem
      nod = elementos(e,n);
      index = zeros(Ndofpornod,1);  
      for i=1:Ndofpornod
        index(i) = (nod-1)*Ndofpornod+i;
      end
      el=(n-1)*Ndofpornod+1:(n)*Ndofpornod;
      elemDof(e,el)= index;
    end
  end
end