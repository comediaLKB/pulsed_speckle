function [TM_cut] = cut_SLM_modes(TM,flag)

%this function aims at cutting the SLM modes outside the pupile defined by
%the backface of the objevtive microscope.
%One has to select the center of the pupile and then one edge.
%This function returns the cutted TM.
%TM should be a 3D matrix, first dimention may be 1, 2nd dimension is SLM
%and 3rd CCD
%Flag is for display or not the cutted part of the modes

n_slm = sqrt(size(TM,2));
TM_norm = zeros(1, n_slm, n_slm);
TM_norm(1,:,:) = reshape(squeeze(mean(abs(TM(1,:,:)).^2,3)), n_slm, n_slm);
%select center and edge
figure(10)
imagesc(squeeze(TM_norm(1,:,:)));
[x_center y_center] = ginput(1);
[x_edge y_edge] = ginput(1);
r = sqrt((x_center-x_edge)^2+(y_center-y_edge)^2);

keep  = [];
TM_cut_plot = TM;
for i=1:size(TM,2)
    idy = mod(i,n_slm);
    idx = (i - mod(i,n_slm))/n_slm;
    a = sqrt((x_center - idx)^2+(y_center - idy)^2);
    if a < r
        keep = [keep i];
    else
        TM_cut_plot(1,i,:) = 0;
    end
end
TM_cut = TM(1,keep,:);

if flag ==1
    TM_norm(1,:,:) = reshape(squeeze(mean(abs(TM_cut_plot(1,:,:)).^2,3)), n_slm, n_slm);
    figure(11)
    imagesc(squeeze(TM_norm(1,:,:)));
end

end