function [gamma] = gamma_evaluation(meta, roi_type)

if roi_type == 'paving'
    gamma = meta.Sub_ROI(1)^2/meta.N^2;
elseif roi_type == 'interlocked'
    gamma = (meta.Sub_ROI(:,1)./N).^2;
else
    gamma = (meta.ROI_tot/N)^2;
end

end