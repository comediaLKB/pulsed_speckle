function [corr] = correlations(A, B)
corr = sum(A.*B)/ sqrt(sum(A.*A)*sum(B.*B));
end