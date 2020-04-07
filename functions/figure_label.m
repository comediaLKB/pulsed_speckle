function [axes_children] = figure_label(axes_parent, number, pos_offset, font_size)

%This code anables to put some labeling text to figures. It takes an axis
%in input and add another one, invisible and with some text, nearby it.

% get axis position
Pos = get(axes_parent,'Position');

% chose default font size from figure height
if nargin < 4
    font_size = 15+2*Pos(4);
end

% text position - default upper edge of the plot
if nargin < 3
    pos_offset = [0 0];
end
text_pos(1) = Pos(1) + pos_offset(1);
text_pos(2) = Pos(2) + Pos(4) + pos_offset(2);
axis_width = 0.3;

% create an invisible axis a bit moved from first one
%axes_children = axes('units', 'centimeters', 'position', [Pos(1)-Pos(3)/4 Pos(2)+Pos(4)+0.25 0.3 0.3]);
axes_children = axes('units', 'centimeters', 'position', [text_pos(1) text_pos(2) axis_width axis_width]);
set(axes_children, 'Visible', 'off')
set(axes_children, 'XTickLabel', [],'YTickLabel', [])

% add some text in it
handle_text = text(0,0,number,'Parent',axes_children,'FontSize', font_size);
set(handle_text,'Interpreter',  'Tex',...
                'FontName',     'Arial',...
                'FontWeight',   'bold',...
                'FontSize',     font_size)

end