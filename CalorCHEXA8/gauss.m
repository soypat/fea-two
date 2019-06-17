function [w, gp, n] = gauss(ng)

dims = length(ng);
n = prod(ng);
w = ones(n,1);
gp = zeros(n,dims);

switch dims

    case 3
        [w1,gp1] = gauss1D(ng(1));
        [w2,gp2] = gauss1D(ng(2));
        [w3,gp3] = gauss1D(ng(3));

        counter = 1;
        for ig1 = 1:ng(1)
            for ig2 = 1:ng(2)
                for ig3 = 1:ng(3)
                    w(counter) = w(counter)*w1(ig1)*w2(ig2)*w3(ig3);
                    gp(counter,:) = [gp1(ig1) gp2(ig2) gp3(ig3)];
                    counter = counter + 1;
                end
            end
        end

    case 2
        [w1,gp1] = gauss1D(ng(1));
        [w2,gp2] = gauss1D(ng(2));

        counter = 1;
        for ig1 = 1:ng(1)
            for ig2 = 1:ng(2)
                w(counter) = w(counter)*w1(ig1)*w2(ig2);
                gp(counter,:) = [gp1(ig1) gp2(ig2)];
                counter = counter + 1;
            end
        end

    case 1
        [w,gp] = gauss1D(n);
        w = w';
        gp = gp';
                
end
        
