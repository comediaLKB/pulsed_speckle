function color = tgtm_paper_colors(color_flag)

if color_flag == 'b'            % blue
    color = [31 119 180]./255;
elseif strcmp(color_flag, 'o')  % orange
    color = [225 127 14]./255;
elseif strcmp(color_flag, 'g')  % green
    color = [44 160 44]./255;
elseif strcmp(color_flag, 'r')  % red
    color = [214 39 40]./255;
elseif strcmp(color_flag, 'y')  % yellow   
    color = [238 178 0]./255;
elseif strcmp(color_flag, 'p')  % purple (4th matlab standard)
    color = [126 47 142]./255;
elseif color_flag == 't'        % teal
    color = [0 150 136]./255;
else
    disp('no standard color!')
    color = [0 0 0];
end

end