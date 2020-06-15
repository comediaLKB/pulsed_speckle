function out_TM = randTM(n_ccd, n_slm, smear_param, form_flag, boundary_flag)

% default parameters
if nargin < 3
    smear_param = 0;    % no grains
    form_flag = 0;      % Gaussian smear
    boundary_flag = 0;  % perodic smear
elseif nargin < 4
    form_flag = 0;      % Gaussian smear
    boundary_flag = 0;  % perodic smear
elseif nargin < 5
    boundary_flag = 0;  % perodic smear
end

out_TM = zeros(n_ccd, n_slm);

if boundary_flag
    roi_vec = (sqrt(n_ccd)/2+1):(1.5*sqrt(n_ccd));
    n_ccd = 4*n_ccd;
end

k = (1:sqrt(n_ccd)) - sqrt(n_ccd)/2 - 1;

rand_TM = randn(n_ccd, n_slm) + 1i * randn(n_ccd, n_slm);

if smear_param
    
    x = (1:sqrt(n_ccd)) - round(sqrt(n_ccd)/2);
    [X,Y] = meshgrid(x,x);
    
    if form_flag == 0 % Gaussian - smear_param is the FWHM grain size
        sigma_grain = smear_param / (2*sqrt(2*log(2)));
        smear_kernal = exp(-(X.^2+Y.^2)/(2*sigma_grain^2));
        smear_kernal = smear_kernal / sqrt(sum(sum(smear_kernal.^2)));
        smear_kernal_k = fft2(smear_kernal);
    elseif form_flag == 1 % Top-hat - smear_param is the maximal k vector
        smear_kernal_k = zeros(sqrt(n_ccd),sqrt(n_ccd));
        for idx_x = 1:sqrt(n_ccd)
            for idx_y = 1:sqrt(n_ccd)
                if sqrt(k(idx_x)^2+k(idx_y)^2) <= smear_param
                    smear_kernal_k(idx_x,idx_y) = 1;
                end
            end
        end
        smear_kernal_k = smear_kernal_k / sum(sum(smear_kernal_k));
        smear_kernal_k = fftshift(smear_kernal_k);
    elseif form_flag == 2 % Goodman - smear_param is the maximal k vector
        % CAREFULL - This is only the sqrt of the Goodman expression which is calculated for the intesity. 
        % Not the proper expression for the field.
        [k_X, k_Y] = meshgrid(k,k);
        k_abs_tilde = sqrt(k_X.^2+k_Y.^2) / smear_param;
        smear_kernal_k = sqrt(real(acos(k_abs_tilde) - k_abs_tilde.*sqrt(1 - k_abs_tilde.^2)));
        smear_kernal_k = smear_kernal_k / sqrt(sum(sum(smear_kernal_k.^2)));
        smear_kernal_k = fftshift(smear_kernal_k);
    end
    
    for idx_slm = 1:n_slm       
        out_line = reshape(rand_TM(:,idx_slm), sqrt(n_ccd), sqrt(n_ccd));
        out_line_smear = fftshift( ifft2( smear_kernal_k .* fft2(out_line) ) );
        
        if boundary_flag
            out_TM(:,idx_slm) = reshape(out_line_smear(roi_vec,roi_vec),n_ccd/4,1);
        else
            out_TM(:,idx_slm) = reshape(out_line_smear,n_ccd,1);
        end
    end
  
else
    
    out_TM(:,:) = rand_TM;
    
end


end