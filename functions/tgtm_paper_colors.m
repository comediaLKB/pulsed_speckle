function color = tgtm_paper_colors(color_flag)

if color_flag == 'b'
    color = [31 119 180]./255;
elseif strcmp(color_flag, 'o')
    color = [225 127 14]./255;
elseif strcmp(color_flag, 'g')
    color = [44 160 44]./255;
elseif strcmp(color_flag, 'r')
    color = [214 39 40]./255;
else
    disp('no standard color!')
    color = [0 0 0];
end

end