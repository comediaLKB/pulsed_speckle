function [Corr] = corr_field(Vector1, Vector2, flag)

%normalization
if flag == 1
    %normalization
    Vector1_norm = Vector1./sqrt(ctranspose(Vector1)*Vector1);
    Vector2_norm = Vector2./sqrt(ctranspose(Vector2)*Vector2);
    %correlation
    Corr = abs(ctranspose(Vector1_norm)*Vector2_norm);
else
%     %normalization
    Vector1_norm = Vector1./sum(Vector1);%(ctranspose(Vector1)*Vector1);
    Vector2_norm = Vector2./sum(Vector2);%(ctranspose(Vector2)*Vector2);
%     %correlation
     Corr = sum((Vector1_norm - mean(Vector1_norm)).*(Vector2_norm - mean(Vector2_norm)));
end
end