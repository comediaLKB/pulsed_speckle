function [focus_vec, focus_point] = focus_center(TM)

n_ccd = size(TM, 1);

focus_point = [sqrt(n_ccd)/2 sqrt(n_ccd)/2];
focus_point_idx = (focus_point(2)-1)*sqrt(n_ccd) + focus_point(1);

focus_vec = exp(1i * angle(TM(focus_point_idx,:)'));

end