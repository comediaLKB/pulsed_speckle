function out = corr_scaled(A,B)
% correlation function scaled by the mean. used in:
% van Beijnum et al. Optics Letters, 36, 373 (2011)
% de Boer et al. PRB, 45, 685 (1992) 

out = mean2( A .* B ) ./ (mean2(A) .* mean2(B)) - 1;


end