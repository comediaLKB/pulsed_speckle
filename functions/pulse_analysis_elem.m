function [images_pulse_ifft_avg, images_pulse_ifft_avg_imag, images_pulse_ifft_avg_imag_field, t, idx] = pulse_analysis_elem(images_all, Ref_all, index_scan, type_idx, meta)

c = 299792458; % speed of light
images_pulse_ifft_avg = zeros(1,size(images_all,2));

    images = double(squeeze(images_all(index_scan,:,:,:)));
    Ref.ref = double(squeeze(Ref_all.ref(index_scan,:,:)));
    Ref.focus = double(squeeze(Ref_all.focus(index_scan,:,:)));
    
    images_pulse = zeros(size(images));
    images_pulse_ifft_avg_imag = zeros(size(images));
    images_pulse_ifft_avg_imag_field = zeros(size(images));
    for t_idx = 1:meta.im
        images_pulse(t_idx,:,:) = (squeeze(images(t_idx,:,:)) - Ref.ref - Ref.focus)./(2*sqrt(Ref.ref .* Ref.focus));
    end
    
    % Axis
    x = linspace(meta.start, meta.stop, meta.im);
    t = 2*(x-x(1))*1e-3 / c;
    t_tot = t(end) - t(1);
    nu = ((1:meta.im) - meta.im/2) * (1/t_tot);

    % Filter
    filter = [310 440];  % in THz
    filter_mask = zeros(size(images_pulse,1),1);
    filter_mask((filter(1) < nu*1e-12) & (filter(2) > nu*1e-12)) = 1;
    
    % Index
    [idx] = reconstruct_idx(type_idx, images_all, index_scan, meta);
%     dist_cam = 3;
% %     if index_scan > size(meta.idx,1)
%         idx(1) = 1;
%         idx(2) = size(images_pulse,2);
%         idx(3) = 1;
%         idx(4) = size(images_pulse,3);
%     else
%         idx = dist_cam*meta.idx(index_scan,:);
%     end
%         
%     if 1; %TMfocus == size(images_all,1)
%         idx(1) = 1;
%         idx(2) = size(images_pulse,2);
%         idx(3) = 1;
%         idx(4) = size(images_pulse,3);
%     else
%         dist_cam = 3;%*meta.dist;
% %         Roi_sub_idx_cam = dist_cam.*Roi_sub_idx;
%         idx = meta.idx(TMfocus,:).*dist_cam;
%     end
    
%     if TMfocus == size(images_all,1)
%         idx(1) = 1;
%         idx(2) = size(images_pulse,2);
%         idx(3) = 1;
%         idx(4) = size(images_pulse,3);
%     else
%          dist_cam = 3*meta.dist;
%          idx(1) = 1+dist_cam*(TMfocus-1);
%          idx(2) = size(images_pulse,2) - dist_cam*(TMfocus-1);
%          idx(3) = 1+dist_cam*(TMfocus-1);
%          idx(4) = size(images_pulse,3) - dist_cam*(TMfocus-1);
% %         sub_roi_num = sqrt(size(images_all,1)-1);
% %         sub_roi_length = size(images_all,3) / sub_roi_num;
% %         
% %         idx(1) = (1 + sub_roi_length*floor((TMfocus-1)/sub_roi_num));
% %         idx(2) = idx(1) + sub_roi_length - 1;
% %         idx(3) = (1 + sub_roi_length*mod((TMfocus-1),sub_roi_num));
% %         idx(4) = idx(3) + sub_roi_length - 1;
%     end
    
    counter = 0;
    for j = idx(1):idx(2) %(Xf-roi_width):(Xf+roi_width)
        for k = idx(3):idx(4) %(Yf-roi_width):(Yf+roi_width)
            
            counter = counter + 1;
            
            images_pulse_fft = fftshift(fft(images_pulse(:,j,k)));
            images_pulse_fft = images_pulse_fft .* filter_mask;
            images_pulse_ifft = ifft(ifftshift(images_pulse_fft));
            
            images_pulse_ifft_avg = images_pulse_ifft_avg + abs(images_pulse_ifft)';
            
             images_pulse_ifft_avg_imag(:,j,k) = abs(images_pulse_ifft);
             images_pulse_ifft_avg_imag_field(:,j,k) = (images_pulse_ifft);
            
        end
    end
    
   images_pulse_ifft_avg = images_pulse_ifft_avg / counter;

end