 function out = corr_mickael(A, B, sub_mean_flag)
 % as used by Mickael after he realized that corr2 gives weird results for
 % the focus. Not normalized.
 
 if nargin < 3
     sub_mean_flag = 0;
 end
 
 if sub_mean_flag
     out = sum(sum( (A - mean(A(:))) .* (B - mean(B(:))) ));
 else
     out = sum(sum( (A) .* conj(B) ));
 end

 end
