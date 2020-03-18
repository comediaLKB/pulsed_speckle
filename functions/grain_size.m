function [size_grain, size_grain_xy, auto_corr_full] = grain_size(image_array, plot_flag, data_flag)
% function to get the grain size of an intesity speckle image or array of images
% BR 09/2019

if nargin < 3
    data_flag = 0;
    if nargin < 2
        plot_flag = 0;
    end
end

size_grain = zeros(size(image_array,1),1);
size_grain_xy = zeros(size(image_array,1),2);
auto_corr_full = zeros(size(image_array,1), size(image_array,2), size(image_array,3));
    
for t_idx = 1:size(image_array,1)

    % get atocorrelation
    image = squeeze(image_array(t_idx,:,:));
    [n, m] = size(image);
    auto_corr = ( fftshift( ifft2(fft2(image) .* conj(fft2(image)) ) ) )./(n*m);

    auto_corr = ( auto_corr - mean(image(:)).^2 ) / mean(image(:)).^2;
    
%     auto_corr = abs(auto_corr);
    
    if data_flag
        auto_corr_full(t_idx,:,:) = auto_corr;
    end
    
%     image2 = image.^2;
%     auto_corr2 = ( fftshift( ifft2(fft2(image2) .* conj(fft2(image2)) ) ) )./(n*m);
%     auto_corr2 = ( auto_corr2 - mean(image2(:)).^2 ) / mean(image2(:)).^2;
    
    if plot_flag
        figure(71)
        clf
        box on
        surf(auto_corr)
    end
    
    % Antonines code
      % normalize
%     auto_corr = auto_corr - min(auto_corr(:));
%     auto_corr = auto_corr / max(auto_corr(:));

%     [value_max, index]=max(auto_corr(:));
%     [x, y]=ind2sub(size(auto_corr),index);
%     xx=numel(auto_corr(x,auto_corr(x,:)>=value_max/2));
%     yy=numel(auto_corr(auto_corr(:,y)>=value_max/2,y));
%     size_grain_a(t_idx) = (xx+yy)/2;
    
    center_x = floor((size(auto_corr,1)-1)/2 +1);
    center_y = floor((size(auto_corr,2)-1)/2 +1);
    x = (-(size(auto_corr,1)-1)/2):((size(auto_corr,1)-1)/2);
    y = (-(size(auto_corr,1)-1)/2):((size(auto_corr,1)-1)/2);
    
    ftype = fittype('c + b*exp(-x.^2 / (2*a^2))','coefficients',{'a','b','c'});
    fo = fitoptions('Method', 'NonlinearLeastSquares','StartPoint',[3, auto_corr(center_x,center_y), 0]);
    gfit_x = fit(x', auto_corr(:,center_y),ftype,fo);
    gfit_y = fit(y', auto_corr(center_x,:)',ftype,fo);

    if plot_flag
        figure(72)
        clf
        box on
        hold on
        plot(x,auto_corr(:,center_y), '.-')
        plot(x,auto_corr(center_x,:), '.-')
%         plot(x,auto_corr2(:,center_y), '.')
        plot(gfit_x)
        plot(gfit_y)
        title(sprintf('idx %d',t_idx))
    end
    
    % FWHM from Gauss sigma
    size_grain_xy(t_idx,1) = 2*sqrt(2*log(2)) * gfit_x.a;
    size_grain_xy(t_idx,2) = 2*sqrt(2*log(2)) * gfit_y.a;
    size_grain(t_idx) = mean(size_grain_xy(t_idx,:)); 
    
end

end