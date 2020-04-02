function [Corr_norm_avg, Delta_r] = Sim_avg_spatial_corr(TM,N_avg,flag)

% takes a TM, apply a random vector on it and compute spatial correlations
% + average according to N_avg value. Flag is for printing averaging
% increase
% This function calls the function spatial_Corr
for i=1:N_avg
    if flag
        i
    end
    Input = randn(size(TM,1),1) + 1i*randn(size(TM,1),1);
    Output = TM'*Input;
    Output = reshape(Output,sqrt(size(TM,2)),sqrt(size(TM,2)));
    [Corr_norm, Delta_r] = Spatial_corr(abs(Output));
    if i == 1
        Corr_norm_avg = zeros(size(Corr_norm));
    end
    Corr_norm_avg = Corr_norm_avg + Corr_norm;
end
Corr_norm_avg = Corr_norm_avg/N_avg;


end