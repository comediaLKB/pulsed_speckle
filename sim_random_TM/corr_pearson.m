function [Corr] = corr_pearson(Vector1, Vector2)

% mean values
Vector1_no_mean = Vector1 - mean(Vector1);
Vector2_no_mean = Vector2 - mean(Vector2);

% correlation
Corr = corr_field(Vector1_no_mean, Vector2_no_mean);

end