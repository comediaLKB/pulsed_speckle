function [Vector_filt] = filtre(Vector, bandwidth, flag)

% fft
fft_corr = fftshift(fft(Vector));

% create the mask filter
filt_mask = zeros(size(fft_corr));
filt_mask(bandwidth)= 1;

%filter
fft_corr_filt = fft_corr.*filt_mask;

if flag == 1
    hold on
    plot(abs(fft_corr))
    plot(filt_mask*max(fft_corr))
end

Vector_filt = ifftshift(ifft(fft_corr_filt));

end