function [out, t, dt] = randTM_dyn(n_ccd, n_slm, n_t, t_max)

% time axis
t = linspace(0, t_max, n_t);
dt = t(2) - t(1);

% gaussian smear function (sigma = 1)
gauss_t = exp(-(t - (t(end/2+1))).^2 / 2);
gauss_t = gauss_t / (sum(gauss_t)*dt);

TM_random = randn(n_t, n_ccd, n_slm) + 1i * randn(n_t, n_ccd, n_slm);

% for idx = 1:length(t)
%     TM_random(idx,:,:) = make_unitary(squeeze(TM_random(idx,:,:)));
% end

out = fftshift(ifft(fft(gauss_t') .* fft(TM_random,[],1),[],1),1) * dt;

% for idx = 1:length(t)
%     out(idx,:,:) = make_unitary(squeeze(out(idx,:,:)));
% end

end