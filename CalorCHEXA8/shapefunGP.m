function [Ns,dNs,dNauxs,NLs] = shapefunGPM(upg,N,dN,dNaux,NL)
npg=size(upg,1);
Ns=cell(npg,1);
dNs = cell(npg,1);
dNauxs = cell(npg,1);
NLs = cell(npg,1);

for ipg =1:npg %Optimiza para problemas grandes
    ksi = upg(ipg,1); eta = upg(ipg,2);zeta = upg(ipg,3);
    
    Ns{ipg} = eval(subs(N));
    dNs{ipg} = eval(subs(dN));
    dNauxs{ipg} = eval(subs(dNaux));
    NLs{ipg} = eval(subs(NL));
end
end