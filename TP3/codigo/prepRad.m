%Preparar radiación
Tespacio = 2.7;
Trad = Tespacio.^4;
boltz = 1.417e-8;
supnod = [1 2 3 4
          1 5 8 4
          4 8 7 3
          3 7 6 2
          2 6 5 1
          5 6 7 8];
k=1/sqrt(3);
wpgs=ones(4,1);
upg2d=[-k -k;-k k;k k;k -k];
upg1 = [-k -k k k]';
upg2 = [-k k k -k]';
upgS = ones(4,1);
npgs=4;
olds =0;
Ns = cell(6,1);
dNauxs = cell(6,1);
for s=1:6
    switch s
        case 1
            surf_upg = [upg1 upg2 -upgS];
        case 2
            surf_upg = [upg1 -upgS upg2];
        case 3
            surf_upg = [upgS upg1  upg2];
        case 4
            surf_upg = [upg1 upgS upg2];
        case 5
            surf_upg = [-upgS upg1  upg2];
        case 6
            surf_upg = [upg1 upg2 upgS];
    end
    [Nm,dNauxm] = shapefunGPM(surf_upg,N,dNaux);
    Ns{s} = Nm;
    dNauxs{s} = dNauxm; %fijate esta optimizacion papu. 2 minutos a 0.3 segundos
    %chau, te kiero elementos finulis. hola metodos numericos II
end