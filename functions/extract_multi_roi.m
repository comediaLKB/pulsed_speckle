function [im, im_ref] = extract_multi_roi(Images_set, meta)

%Images_set is of the size n,m,m; m being the full ROI and n the number of
%realizations = meta.N_roi

% extract values from Roi
im =  zeros(meta.N_roi,size(Images_set,2)/sqrt(meta.N_roi),size(Images_set,3)/sqrt(meta.N_roi));
im_ref = zeros((meta.N_roi-1)*meta.N_roi,size(Images_set,2)/sqrt(meta.N_roi),size(Images_set,3)/sqrt(meta.N_roi));
curr_idx = 0;
for i=1:meta.N_roi
    A = squeeze(Images_set(i,:,:));
    idx = [1+meta.bin_val*(meta.idx(i,1)-1) meta.bin_val*meta.idx(i,2) 1+meta.bin_val*(meta.idx(i,3)-1) meta.bin_val*meta.idx(i,4)];
    im(i,:,:) = A(idx(1):idx(2),idx(3):idx(4));
    for j=1:meta.N_roi
        if i == j
        else
            curr_idx = curr_idx +1;
            idx = [1+meta.bin_val*(meta.idx(j,1)-1) meta.bin_val*meta.idx(j,2) 1+meta.bin_val*(meta.idx(j,3)-1) meta.bin_val*meta.idx(j,4)];
            im_ref(curr_idx,:,:) = A(idx(1):idx(2),idx(3):idx(4));
        end
    end
end
end