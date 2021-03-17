function [TM_filt, TM_pic_fft] = filter_TM(TM, filter_width, filter_flag, ring_width)

if nargin < 4
    ring_width = 1;
end

[n_ccd, n_slm] = size(TM);

k = (1:sqrt(n_ccd)) - sqrt(n_ccd)/2 - 1;
[KX,KY] = meshgrid(k,k);

TM_filt = zeros(size(TM));
TM_pic_fft = zeros(sqrt(n_ccd),sqrt(n_ccd));
for idx_slm = 1:n_slm

    TM_pic = reshape(TM(:,idx_slm), sqrt(n_ccd), sqrt(n_ccd));
    
    if filter_flag == 1 % 1D Gaussian smear
        smear_kernal_k = exp(-KX.^2 / (2*filter_width^2));
        smear_kernal_k = smear_kernal_k / sum(sum(smear_kernal_k));
        smear_kernal_k = fftshift(smear_kernal_k);
        TM_pic = fftshift( ifft2( smear_kernal_k .* fft2(TM_pic) ) );        
    elseif filter_flag == 2 % ring
        smear_kernal_k = zeros(size(TM_pic));
        k_min = filter_width - ring_width;
        k_max = filter_width;
        smear_kernal_k((sqrt(KX.^2 + KY.^2) >= k_min) & (sqrt(KX.^2 + KY.^2) <= k_max)) = 1;
        smear_kernal_k = smear_kernal_k / sum(sum(smear_kernal_k));
        smear_kernal_k = fftshift(smear_kernal_k);
        TM_pic = ( ifft2( smear_kernal_k .* fft2(TM_pic) ) );
    elseif filter_flag == 3 % phase vortex
        shift_angle = -atan(KY./KX)+pi/2;
        shift_angle(:,1:(sqrt(n_ccd)/2)) = shift_angle(:,1:(sqrt(n_ccd)/2)) + pi;
        smear_kernal_k = exp(1i*shift_angle) / n_ccd;
        smear_kernal_k((sqrt(n_ccd)/2+1),(sqrt(n_ccd)/2+1)) = 0;
        smear_kernal_k = fftshift(smear_kernal_k);
        TM_pic = fftshift( ifft2( smear_kernal_k .* fft2(TM_pic) ) );
    elseif filter_flag == 4 % single k-pixel   
        smear_kernal_k = zeros(size(TM_pic));
        smear_kernal_k(k==filter_width,k==filter_width) = 1;
        smear_kernal_k = fftshift(smear_kernal_k);
        TM_pic = ( ifft2( smear_kernal_k .* fft2(TM_pic) ) );
    end

    TM_pic_fft = TM_pic_fft + abs(fftshift(fft2(TM_pic)))/n_slm;

    TM_filt(:,idx_slm) = reshape(TM_pic, n_ccd, 1);

end

end