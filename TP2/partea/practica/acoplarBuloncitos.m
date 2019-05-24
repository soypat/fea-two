function [Kg fixity] = acoplarBuloncitos(azules,K,kb)
% LEGACY. No usar
n2d6=@(n) [n.*6-5 n.*6-4 n.*6-3 n.*6-2 n.*6-1 n.*6];
n2d3 = @(n) [n.*3-2 n.*3-1 n.*3];

Nazul = length(azules);
extranod = Nazul * 3;
extradof = extranod*3;
Ndofpornodo = 6;
Nnodi = size(K)/6;
elemDofEnd = n2d6(azules);
dof = size(K,1);
Kg = zeros(dof+extradof,dof+extradof);
Kg(1:dof,1:dof)=K;

giroDofs = [0 0 0 1 1 1]; % No quiero otorgar rigidez a los giros.
fixity = false(dof+extradof,1);
for kazul=1:length(azules)
    storeToAzul = n2d6(azules(kazul));
    storeToAzul = storeToAzul(~giroDofs);
    Kbx = vigastiffness(kb(1),0,0,1,0,0,0,1); 
    Kby = vigastiffness(kb(2),0,0,1,0,0,0,1);
    Kbz = vigastiffness(kb(3),0,0,1,0,0,0,1);
    Kbx = vigorotar(Kbx,[0 0 0],[1 0 0],[0 0 1]);
    Kby = vigorotar(Kby,[0 0 0],[0 1 0],[0 0 1]);
    Kbz = vigorotar(Kbz,[0 0 0],[0 0 1],[1 0 0]);
    Kbx=Kbx(~[giroDofs giroDofs],~[giroDofs giroDofs]);
    Kby=Kby(~[giroDofs giroDofs],~[giroDofs giroDofs]);
    Kbz=Kbz(~[giroDofs giroDofs],~[giroDofs giroDofs]);
    
    for rl = 0:2
        storeToCB = dof+n2d3(kazul*3-2+rl);
        storeTo = [storeToAzul storeToCB];
        switch rl
            case 0
                Kg(storeTo,storeTo) = Kg(storeTo,storeTo)+Kbx;
                fixity(storeTo) = fixity(storeTo) | [0 0 0 1 0 0]';
            case 1
                Kg(storeTo,storeTo) = Kg(storeTo,storeTo)+Kby;
                fixity(storeTo) = fixity(storeTo) | [0 0 0 0 1 0]';
            case 2
                Kg(storeTo,storeTo) = Kg(storeTo,storeTo)+Kbz;
                fixity(storeTo) = fixity(storeTo) | [0 0 0 0 0 1]';
        end
    end
end

end

