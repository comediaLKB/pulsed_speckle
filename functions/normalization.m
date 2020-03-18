function [Image_norm] = normalization(Image, Norm, view)

% this works assuming size(images) = columns lines

Image_norm = zeros(size(Image));

 if isequal(Norm ,'lines')
     Image_norm= Image ./ sqrt(sum(Image .* conj(Image)));

 elseif isequal(Norm ,'columns')
     Image_temp = transpose(Image);
     Image_temp = Image_temp ./ sqrt(sum(Image_temp .* conj(Image_temp)));
     Image_norm = transpose(Image_temp);
     
 elseif isequal(Norm ,'lines_columns')
     Image_temp = Image ./ sqrt(sum(Image .* conj(Image)));
     Image_temp = Image_temp' ./ sqrt(sum(Image_temp' .* conj(Image_temp')));
     Image_norm = Image_temp';
     
 else
     disp('No Normalization defined here')
 end
 
 if view == 1
     imagesc(abs(Image_norm));
 end

end
