 function [corr] = correlations_mat(A, B)
% corr = sum(sum(A.*B))/ sqrt(sum(sum(A.*A))*sum(sum(B.*B)));
% corr = 1 - mean2(abs(A-B))/sqrt(mean2(B)*mean2(A));
corr = sum(sum((A-mean2(A)).*conj(B-mean2(B))));
%corr = sum(sum(A.*conj(B)));
end