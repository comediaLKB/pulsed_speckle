function h = imagecf(varargin)
%IMAGECF    Display complex field 
%   IMAGECF(Z) 
%   IMAGECF(Z,Zmax) 
%   IMAGECF(X+iY,Z) 
%   IMAGECF(X+iY,Z,Zmax) 
%   IMAGECF(ax,___) 
%   IMAGECF(___, option) 
%       color of zero.
%          'w' white
%          'k' black
%   Example>
%       [x,y] = meshgrid(-2:.1:2);
%       z = x + 1i*y;
%       imagecf(x+1i*y, z.^2+1,1,'k');
% 
%   The code is firstly written by Kanghyun Chu. (kanghyunchu@kaist.ac.kr)
%   2018-Apr-10 Version 1.0

%%
[cax,args,nargs] = axescheck(varargin{:});
%   cax : obtained axes object
%   args: arguments except axes handle
%   nargs = numel(args)
%% close previous colormap
fo = findobj('type','figure');
for ii = 1:numel(fo)
    if strcmp(fo(ii).Name,'Color scale map')
        try
            close(fo(ii));
        end
        break
    end
end
clear fo ii
%%
switch nargs
    case 1
%   IMAGECF(Z) 
        if isnumeric(args{1})
            Z = args{1};
            [X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
            zc = [];
            Z_max = [];
        else 
            error(message('MATLAB:numericArgRequired'));
        end
        
    case 2
        if ischar(args{2})
% Z, 'zc'
            Z = args{1};
            [X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
            zc = args{2};
            Z_max = [];
        elseif numel(args{2}) == 1
% Z, Zmax
            Z = args{1};
            [X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
            zc = [];
            Z_max = args{2};
        else
% X+iY, Z
            Z = args{2};
            X = real(args{1});  Y = imag(args{1});
            zc = [];
            Z_max = [];
            if ~isequal(size(Z), size(X))
                error(message('MATLAB:arrayDimensionsMustMatch')); 
            end
        end
        
    case 3
        if ischar(args{3})
% Z, Zmax, 'zc'
            if numel(args{2}) == 1
                Z = args{1};
                [X,Y]=meshgrid(1:size(Z,2),1:size(Z,1));
                zc = args{3};
                Z_max = args{2};
            else
% X+iY, Z, 'zc'
                Z = args{2};
                X = real(args{1});  Y = imag(args{1});
                zc = args{3};
                Z_max = [];
                if ~isequal(size(Z), size(X))
                    error(message('MATLAB:arrayDimensionsMustMatch')); 
                end       
            end
            
        elseif numel(args{3}) == 1
% X+iY, Z, Zmax
            Z = args{2};
            X = real(args{1});  Y = imag(args{1});
            zc = [];
            Z_max = args{3};
            if ~isequal(size(Z), size(X))
                error(message('MATLAB:arrayDimensionsMustMatch')); 
            end
        else
            error(message('MATLAB:badargs'));
        end
        
    case 4
% X+iY, Z, Zmax, 'zc'
        Z = args{2};
        X = real(args{1});  Y = imag(args{1});
        zc = args{4};
        Z_max = args{3};
        
    otherwise
        error(message('MATLAB:narginchk:tooManyInputs'));
end
%% default values
if isempty(zc)
    zc = 'w';
end
if isempty(Z_max)
    Z_max = max( abs(Z(:)));
end

%% zero color - default: white
if zc == 'k'
        zclr = [0 0 0];
elseif zc == 'w'
        zclr = [1 1 1];
else
        zclr = [1 1 1];
        warning(message('MATLAB:badopt'))
end

%% colorspace define 
% clr = hsv(128); % hsv color
clr = 1/2 + sqrt(6)/4*(kron([2;-1;-1]/sqrt(6), cos( linspace(0,2*pi*127/128,128) ))...
                      +kron([0; 1;-1]/sqrt(2), sin( linspace(0,2*pi*127/128,128) ))); % intensity normalized rgb color ring
clr = clr'; 
cmap = zeros(128,128,3); % Order: amp, phase, rgb space 
cmap(end,:,:) = reshape(clr,1,[],3);
cmap(1,:,:)   = repmat( reshape(zclr,1,1,3), 1,128);

for i = 1:128
    for k=1:3
        cmap(:,i,k) = linspace(cmap(1,i,k),cmap(end,i,k), 128)';
    end
end
cmap = reshape(cmap,128^2,3);

%% corresponding rgb index
inZa = floor( abs(Z)/Z_max * 127)+1;
    inZa(inZa > 128) = 128; % out-of-range saturated
inZph= mod(round(angle(Z)/(2*pi)* 128),128)+1;
% amp   1~128  
% phase 1~128  1: -0.5d¥è/2 ~ +0.5d¥è, 
%              2: +0.5d¥è/2 ~ +1.5d¥è ...
inZ = sub2ind([128,128], inZa, inZph);
%% create plot
if isempty(cax)
    f=gcf; 
    newplot(f.CurrentAxes);
else 
    cax=newplot(cax);
    cax.Parent.CurrentAxes = cax;
end
% h = warp(X,Y,zeros(size(X)),inZ, cmap);
h = warp(X,Y,inZa,inZ, cmap); % surf Z height = amplitude
    h.Parent.View = [0 90];
    h.Parent.YDir ='normal';
pos  = h.Parent.Parent.Position;

%% colormap display - overwrite openned colormap
% f_ref=figure('Name','Color scale map','Position',[pos(1)+pos(3)+5, pos(2)+pos(4)-320, 320, 320],'Color','w'); % colormap display
% ax_ref = axes(f_ref,'Position',[0.05 0.01 0.90 0.90]);
% [r,t] = meshgrid(1:128, mod(-64:64,128)+1);
% Z_ref = sub2ind([128,128], r, t);
% 
% warp(r/128.*cos(t*2*pi/128), r/128.*sin(t*2*pi/128),r/128, Z_ref,cmap);
% 
% set(ax_ref,   'View',[0 90], 'YDir','normal', 'DataAspectRatio', [1 1 1],'XLim',[-1.5 1.5],'YLim',[-1.5 1.5], ...
%               'XTick', [-1 0 1], 'YTick', [-1 0 1],'YTickLabel',{'' ,'' ,''}, ...
%               'XAxisLocation','origin','YAxisLocation','origin','Layer','top');
% ax_ref.XTickLabel = {'-|Z_{max}|           ', '     0', '           +|Z_{max}|'} ;
% ax_ref.FontName ='Times new roman';
% ax_ref.FontSize = 12;
% xlabel('Real','FontName','Times new roman','FontSize',14);
% ylabel('Imag','FontName','Times new roman','FontSize',14);
% title('Color scale map');
% text(-1.5, 1.1, sprintf('|Z_{max}| = %0.5G',Z_max),'FontName','Times new roman','FontSize',12);
end
