function [grain_size] = grain_size_from_TM(TM,meta)

% This function returns the grain size from computing the participation
% number from the Transmission Matrix

[U S V] = svd(TM);
Sing_values = diag(S);
Participation_number = sum(Sing_values).^2/sum(Sing_values.^2);
Participation_number_norm = Participation_number/min(size(TM,1),size(TM,2));
grain_size = meta.bin_val/Participation_number_norm;
end