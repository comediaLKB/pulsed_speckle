function output = k_dist(image)

output = fftshift( abs(fft2(image)).^2 ); 

end