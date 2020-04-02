function [Corr_norm, Delta_r] = Spatial_corr(Image)

% This function computes spatial correlations. From an image it computes
% the correlations between pixels of a fixed distance (dr is here the euclidian distance) 

Corr = [];
Corr2 = [];
Delta_r = [];
Counter = [];

for i0=1:size(Image,1)
    for j0=1:size(Image,2)
        temp = Image(i0,j0).*Image;
        for i=1:size(Image,1)
            for j=1:size(Image,2)
                dr = (i0-i)^2 + (j0-j)^2;% sqrt taken later to prevent computationnal errors
                [Lia, Locb] = ismember(dr,Delta_r);
                if not(Lia)
                    Corr = [Corr 0];
                    Corr2 = [Corr2 0];
                    Counter = [Counter 0];
                    Delta_r = [Delta_r dr];
                    Locb = size(Delta_r,2);
                end
                Corr(Locb) = Corr(Locb) + temp(i,j);
                Corr2(Locb) = Corr2(Locb) + Image(i0,j0);
                Counter(Locb) = Counter(Locb) + 1;
            end
        end
    end
end

% Sorting and Normalisation
[Delta_r,idx] = sort(Delta_r);
Delta_r = sqrt(Delta_r);
Corr = Corr(idx);
Corr2 = Corr2(idx);
Counter = Counter(idx);

Corr_norm = Corr./(Counter);
Corr2_norm = Corr2./(Counter);
Corr_norm = Corr_norm./Corr2_norm.^2 - 1;

end