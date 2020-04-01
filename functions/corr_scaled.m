function out = corr_scaled(A,B)
% correlation function scaled by the mean. used in ﻿van Beijnum et al. Optics Letters, 36, 373 (2011) 

out = sum(sum( (A / mean2(A)) .* (B / mean2(B)) )) / (size(A,1)*size(A,2)) - 1;

end