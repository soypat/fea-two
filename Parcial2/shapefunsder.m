function [dN, dNinc] = shapefunsder(pointArray,eleType)


ngauss = size(pointArray,1);

switch eleType

    case {'Q9', 'AHMAD9'}
        dN = zeros(2,9,ngauss);
        dNinc = zeros(2,2,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            dN(:,:,igauss) = [ % derivadas respecto de ksi
                (eta*(2*ksi - 1)*(eta - 1))/4,      (eta*(2*ksi + 1)*(eta - 1))/4,       (eta*(2*ksi + 1)*(eta + 1))/4,...
                (eta*(2*ksi - 1)*(eta + 1))/4,                 -eta*ksi*(eta - 1),        -((eta^2 - 1)*(2*ksi + 1))/2,...
                           -eta*ksi*(eta + 1),       -((eta^2 - 1)*(2*ksi - 1))/2,           2*ksi*(eta^2 - 1)
                % derivadas respecto de eta
                  (ksi*(2*eta - 1)*(ksi - 1))/4,    (ksi*(2*eta - 1)*(ksi + 1))/4,       (ksi*(2*eta + 1)*(ksi + 1))/4,...
                  (ksi*(2*eta + 1)*(ksi - 1))/4,     -((2*eta - 1)*(ksi^2 - 1))/2,                -eta*ksi*(ksi + 1),...
                   -((2*eta + 1)*(ksi^2 - 1))/2,               -eta*ksi*(ksi - 1),           2*eta*(ksi^2 - 1) ];
        end
    case {'Q9inc'}
        dN = zeros(2,9,ngauss);
        dNinc = zeros(2,2,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);
            dNinc(:,:,igauss) = [%Derivadas respecto de ksi
                                 3*ksi^2-1      0        
                                 %Derivadas respecto de ksi
                                 0         3*eta^2-1];
            dN(:,:,igauss) = [ % derivadas respecto de ksi
                (eta*(2*ksi - 1)*(eta - 1))/4,      (eta*(2*ksi + 1)*(eta - 1))/4,       (eta*(2*ksi + 1)*(eta + 1))/4,...
                (eta*(2*ksi - 1)*(eta + 1))/4,                 -eta*ksi*(eta - 1),        -((eta^2 - 1)*(2*ksi + 1))/2,...
                           -eta*ksi*(eta + 1),       -((eta^2 - 1)*(2*ksi - 1))/2,           2*ksi*(eta^2 - 1)
                % derivadas respecto de eta
                  (ksi*(2*eta - 1)*(ksi - 1))/4,    (ksi*(2*eta - 1)*(ksi + 1))/4,       (ksi*(2*eta + 1)*(ksi + 1))/4,...
                  (ksi*(2*eta + 1)*(ksi - 1))/4,     -((2*eta - 1)*(ksi^2 - 1))/2,                -eta*ksi*(ksi + 1),...
                   -((2*eta + 1)*(ksi^2 - 1))/2,               -eta*ksi*(ksi - 1),           2*eta*(ksi^2 - 1) ];
        end

    case {'Q8', 'AHMAD8'}
        dN = zeros(2,8,ngauss);
        dNinc = zeros(2,2,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);
            
            dN(:,:,igauss) = [  % derivadas respecto de ksi
                -0.25*(-1+eta)*(eta+2*ksi),  -0.25*(-1+eta)*(-eta+2*ksi),    0.25*(1+eta)*(eta+2*ksi),   0.25*(1+eta)*(-eta+2*ksi),...
                              ksi*(-1+eta),        -0.5*(-1+eta)*(1+eta),                -ksi*(1+eta),        0.5*(-1+eta)*(1+eta)
                % derivadas respecto de eta
                -0.25*(-1+ksi)*(ksi+2*eta),   -0.25*(1+ksi)*(ksi-2*eta),    0.25*(1+ksi)*(ksi+2*eta),   0.25*(-1+ksi)*(ksi-2*eta),...
                      0.5*(-1+ksi)*(1+ksi),                -(1+ksi)*eta,       -0.5*(-1+ksi)*(1+ksi),               (-1+ksi)*eta ];
        end
    case 'Q8inc'
        dN = zeros(2,8,ngauss);
        dNinc = zeros(2,2,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);
            dNinc(:,:,igauss) = [%Derivadas respecto de ksi
                                 3*ksi^2-1      0        
                                 %Derivadas respecto de ksi
                                 0         3*eta^2-1];
            dN(:,:,igauss) = [  % derivadas respecto de ksi
                -0.25*(-1+eta)*(eta+2*ksi),  -0.25*(-1+eta)*(-eta+2*ksi),    0.25*(1+eta)*(eta+2*ksi),   0.25*(1+eta)*(-eta+2*ksi),...
                              ksi*(-1+eta),        -0.5*(-1+eta)*(1+eta),                -ksi*(1+eta),        0.5*(-1+eta)*(1+eta)
                % derivadas respecto de eta
                -0.25*(-1+ksi)*(ksi+2*eta),   -0.25*(1+ksi)*(ksi-2*eta),    0.25*(1+ksi)*(ksi+2*eta),   0.25*(-1+ksi)*(ksi-2*eta),...
                      0.5*(-1+ksi)*(1+ksi),                -(1+ksi)*eta,       -0.5*(-1+ksi)*(1+ksi),               (-1+ksi)*eta ];
        end
    case {'Q4', 'AHMAD4'}
        dN = zeros(2,4,ngauss);
        dNinc = zeros(2,2,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            dN(:,:,igauss) = [  % derivadas respecto de ksi
                -0.25*(1 - eta),  0.25*(1 - eta), 0.25*(1 + eta), -0.25*(1 + eta)
                % derivadas respecto de eta
                -0.25*(1 - ksi), -0.25*(1 + ksi), 0.25*(1 + ksi),  0.25*(1 - ksi) ];
        end

end

