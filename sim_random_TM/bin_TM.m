function [TM_bin, n_ccd_bin] = bin_TM(TM,bin)

[n_ccd, n_slm] = size(TM);
n_ccd_bin = n_ccd / bin^2;

TM_bin = zeros(n_ccd_bin, n_slm);
for idx_slm = 1:n_slm

    TM_pic = reshape(TM(:,idx_slm), sqrt(n_ccd), sqrt(n_ccd));
    
    TM_pic_bin = zeros(sqrt(n_ccd_bin), sqrt(n_ccd_bin));
    for idx_x = 1:sqrt(n_ccd_bin)
        for idx_y = 1:sqrt(n_ccd_bin)
            
            idx_x_range = ((idx_x-1)*bin+1):(idx_x*bin);
            idx_y_range = ((idx_y-1)*bin+1):(idx_y*bin);
            
            TM_pic_bin(idx_x,idx_y) = mean(mean(TM_pic(idx_x_range,idx_y_range)));
            
        end
    end
    
    TM_bin(:,idx_slm) = reshape(TM_pic_bin, n_ccd_bin, 1);

end

end