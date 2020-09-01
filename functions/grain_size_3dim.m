function [grain_size_x, grain_size_y, grain_size_z] = grain_size_3dim(image_array, plot_flag, data_flag)

% this function computes the 3d speckle grain size assuming that 2nd
% dimention is x then y and then z. The first dimension of Image should be
% the averaging.

%initialisation
grain_size_x = 0;
grain_size_y = 0;
grain_size_z = 0;

for t_idx = 1:size(image_array,1)
    t_idx
    
    % speckle grain size in (y,z) plane
    for x_0 = 1:size(image_array,2)
        [size_grain_temp, size_grain_xy_temp, auto_corr_full] = grain_size(reshape(abs(squeeze(image_array(t_idx,x_0,:,:))),1,size(image_array,3),size(image_array,4)), plot_flag, data_flag);
        grain_size_y = grain_size_y + size_grain_xy_temp(1,1);
        grain_size_z = grain_size_z + size_grain_xy_temp(1,2);
    end
    for y_0 = 1:size(image_array,3)
        [size_grain_temp, size_grain_xy_temp, auto_corr_full] = grain_size(reshape(abs(squeeze(image_array(t_idx,:,y_0,:))),1,size(image_array,2),size(image_array,4)), plot_flag, data_flag);
        grain_size_x = grain_size_x + size_grain_xy_temp(1,1);
        grain_size_z = grain_size_z + size_grain_xy_temp(1,2);
    end
    for z_0 = 1:size(image_array,4)
        [size_grain_temp, size_grain_xy_temp, auto_corr_full] = grain_size(reshape(abs(squeeze(image_array(t_idx,:,:,z_0))),1,size(image_array,2),size(image_array,3)), plot_flag, data_flag);
        grain_size_x = grain_size_x + size_grain_xy_temp(1,1);
        grain_size_y = grain_size_y + size_grain_xy_temp(1,2);
    end
end
grain_size_x = grain_size_x/(size(image_array,1)*(size(image_array,3)+size(image_array,4)));
grain_size_y = grain_size_y/(size(image_array,1)*(size(image_array,2)+size(image_array,4)));
grain_size_z = grain_size_z/(size(image_array,1)*(size(image_array,2)+size(image_array,3)));

end