function [M] = vortex_mask(Nv, nbLin, nbCol)
    k = [];
    theta = [];
    for j=1:Nv
        k = [k (2*j-1)/Nv - 1];
        theta = [theta 0];
    end
    M = zeros(nbLin, nbCol);
    for a=1:nbLin
        for b=1:nbCol
            z = (2*b-1)/nbCol - 1 + 1i * (2*a-1)/nbLin - 1i;
            if abs(z) <= 1
                v = 1;
                for j=1:Nv
                    v = v * (z - k(j) * exp(1i * theta(j)));
                end
                M(a,b) = exp(1i * angle(v));
            end
        end
    end 
end