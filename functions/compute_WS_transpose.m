function [Q] = compute_WS_transpose(Ecan_tot,delta_lam,lambda_0,N)

nm = 1E-9;
ps = 1E-12;
c = 299792458;
Q = zeros(N^2, N^2);

% Q(:,:) = 1i*squeeze(Ecan_tot(2,:,:))*(squeeze(Ecan_tot(3,:,:))'-squeeze(Ecan_tot(1,:,:))')/(2*delta_lam*nm)*((lambda_0*nm)^2/(2*pi*c));
% 
Q(:,:) = -1i*squeeze(Ecan_tot(2,:,:))*(squeeze(Ecan_tot(3,:,:))'-squeeze(Ecan_tot(1,:,:))')/(2*delta_lam*ps);

end