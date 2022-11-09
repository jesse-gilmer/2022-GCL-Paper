function H=colored_line_2D(x_coordinates,y_coordinates,color,linestyle,more_varargin)
%---------------------------------------------------------------------------------------------
% H=colored_line_2D(x_coordinates,y_coordinates,color,linestyle,more_varargin)
% This routine draws a line (2d) in a certain 'color' and a certain
% 'linestyle' (which have to be passed as strings, e.g. 'red' (or 'r') and '-',
% or 'green' (or 'g') and ':', etc...)
% More input argument can be passed through 'more_varargin' that must be a cell
% array of valid line options (e.g. {'LineWidth',2})
%
% Inside it uses the basic 'line' function of Matlab (which always gives a
% blue solid line), hence the 'x_coordinates' (of initial and final point of the line),
%'y_coordinates', 'z_coordinates' must be passed just like they would be
%passed to the standard 'line' function
%
%   Usage Examples:
%
%   colored_line([0,3],[0,2],[0,1],'r',':'); 
%   
%   draws a red dotted line between points (0;0;0) and (3;2;1)
%
%   If we want to specify a weight of 2, the above line must be modified to:
%
%   colored_line_2d([0,3],[0,2],'r',':',{'LineWidth',2});
%
%   Matteo Mischiati
%   September 2013
%---------------------------------------------------------------------------------------------

if (nargin <3)
 color='b';
 linestyle='-';
elseif (nargin==3)
    linestyle='-';
end

H=line(x_coordinates,y_coordinates);
set(H,'Color',color,'LineStyle',linestyle);

if nargin > 4
    for ii = 1:2:size(more_varargin, 2)
        set(H, more_varargin{ii}, more_varargin{ii+1});
    end
end