function out = corr_auto(image)

    [n, m] = size(image);
    %out = ( fftshift( ifft2(fft2(image) .* conj(fft2(image)) ) ) )./(n*m);
    out = ( fftshift( ifft2(abs(fft2(image)) ) ) )./(n*m);

    out = ( out - mean(image(:)).^2 ) / mean(image(:)).^2;
    
end