function [Corr] = corr_field(Vector1, Vector2)
 
% normalization
Vector1_norm = Vector1 ./ sqrt( Vector1' * Vector1 );
Vector2_norm = Vector2 ./ sqrt( Vector2' * Vector2 );

% correlation
Corr = Vector1_norm' * Vector2_norm;

end