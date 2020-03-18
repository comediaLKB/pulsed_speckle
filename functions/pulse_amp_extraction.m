function [pulse_amp_avg, pulse_field] = pulse_amp_extraction(images_pulse, t, lambda, filter_width, field_flag)
% function to filter the laser frequency and extract the pulse amplitude and field from a stage trace
% BR 01/2020
%
%   images_pulse.......images for different delay times
%   t..................time axis (s)
%   lambda.............central wavelength of the pulse (m)
%   filter_width.......frequency filter width around pulse center freq (Hz) 
%   field_flag.........should the whole field be saved ? (default = 1)

if nargin < 5
    field_flag = 1;
end

c = 299792458; % speed of light

nr_imag = size(images_pulse,1);
nr_px = size(images_pulse,2)*size(images_pulse,3);

% obtain frequecy axis
t_tot = t(end) - t(1);
nu = ((1:nr_imag) - nr_imag/2) * (1/t_tot);
nu_center = c/lambda; % pulse center frequencz

% define filter mask
filter_mask = zeros(size(images_pulse,1),1);
filter_mask(((nu_center - filter_width) < nu) & ((nu_center + filter_width) > nu)) = 1;
% here take the positive pic

% fft filtering
pulse_amp_avg = zeros(1,size(images_pulse,1));
pulse_field = zeros(size(images_pulse));
for j = 1:size(images_pulse,2)
    for k = 1:size(images_pulse,3)
                
        images_pulse_fft = fftshift(fft(images_pulse(:,j,k)));
        images_pulse_fft = images_pulse_fft .* filter_mask;
        images_pulse_ifft = ifft(ifftshift(images_pulse_fft));
        
        pulse_amp_avg = pulse_amp_avg + abs(images_pulse_ifft)' / nr_px;
        
        if field_flag
            pulse_field(:,j,k) = conj(images_pulse_ifft);
            %here conj is taken because the positive pic is filterred in
            %order to get the same phase as if a phase stepping was
            %performed
        end
        
    end
end
        
end