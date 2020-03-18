function P = exp_dist_modify(I, I_0, g)

x = I / I_0;

P = exp(-x) .* ( 1 + (1/(3*g)) * (x.^2 - 4*x + 2) );
%P = exp(-x);

end