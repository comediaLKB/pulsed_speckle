function k_dist_avg = k_dist_avg(images, sub_mean_flag)

k_dist_indiv = zeros(size(images));
for idx = 1:size(images,1)
    if sub_mean_flag
        images(idx,:,:) = images(idx,:,:) - mean(mean(images(idx,:,:)));
    end
    
    k_dist_indiv(idx,:,:) = k_dist(squeeze(images(idx,:,:)));
    
end

k_dist_avg = squeeze(mean(k_dist_indiv,1));

end