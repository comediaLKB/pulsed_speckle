function out = svd_vector(TM,idx)

[~, ~, V] = svd(TM);
out = V(:,idx);

end