function [axes_children] = figure_label(axes_parent,number)

%This code anables to put some labeling text to figures. It takes an axis
%in input and add another one, invisible and with some text, nearby it.

% get axis position
Pos = get(axes_parent,'Position');

% create an invisible axis a bit moved from first one
axes_children = axes('units', 'centimeters', 'position', [Pos(1)-Pos(3)/4 Pos(2)+Pos(4)+0.25 1 1])
set(axes_children, 'Visible', 'off')

% add some text in it
text(0.25,0.25,number,'Parent',axes_children,'FontWeight','bold','FontSize',15+2*Pos(4));
end