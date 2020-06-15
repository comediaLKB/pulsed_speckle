function out = rand_vector(n)

out = randn(n, 1) + 1i * randn(n, 1);

out = out / sqrt( sum( abs(out).^2 ) );

end