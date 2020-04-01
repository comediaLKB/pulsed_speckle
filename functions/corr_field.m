function [Corr] = corr_field(Vector1, Vector2, norm_flag)

if nargin < 3
    norm_flag = 1;
end


if norm_flag == 1 % field correlation as used (?) in Ambichel ﻿PRX 7, 041053 (2017)
    
    % normalization
    Vector1_norm = Vector1 ./ sqrt( Vector1' * Vector1 );
    Vector2_norm = Vector2 ./ sqrt( Vector2' * Vector2 );

    % correlation
    Corr = real(Vector1_norm' * Vector2_norm);
        
elseif norm_flag == 2 % phase only correlation 
    
    % normalization to phase only
    Vector1_norm = Vector1 ./ abs( Vector1 );
    Vector2_norm = Vector2 ./ abs( Vector2 );
    
    % correlation
    Corr = real(Vector1_norm' * Vector2_norm) / length(Vector1);
    
else
    
    % normalization
    Vector1_norm = Vector1 ./ sum(Vector1);
    Vector2_norm = Vector2 ./ sum(Vector2);
    
    % correlation
    Corr = sum((Vector1_norm - mean(Vector1_norm)) .* (Vector2_norm - mean(Vector2_norm)));
    
end

end