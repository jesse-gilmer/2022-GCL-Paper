function varargout = raster(st, varargin)

%RASTER   Spike raster plot.
%
%   RASTER(ST) given a cell array of spike times ST generates a raster plot.
%
%   Optional args:
%
%   'tlim'     Time limits
%   'clr'      Color (or colormap if matrix)
%   'order'    Row order in which cells should be displayed
%   'byrate'   Order cells by number of spikes
%   'linewidth' Set linewidth value
%
%   original by Ming Gu, changes by CMW

persistent clr;
persistent perm;

% defaults
th = .75; % tick height 
spc = 1;
linewidth = 1;
doblack = 0;
% default MATLAB color order
%{
clr = [ ...
         0         0    1.0000;  ...
         0    0.5000         0;  ...
    1.0000         0         0;  ...
         0    0.7500    0.7500;  ...
    0.7500         0    0.7500;  ...
    0.7500    0.7500         0;  ...
    0.2500    0.2500    0.2500   ...
    ];
%}

clr = [...
    0         0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    ];
    
byrate = 0;
if isnumeric(st), st = {st}; end
perm = 1:length(st);
tlim = [];
use_3d  = 0;

v = varargin;

assignopt( {'color','clr'}, 'clr', v );
assignbopt( 'black', 'doblack', 'v' );
assignopt( 'linewidth', 'linewidth', v );
assignopt( 'order', 'perm', v );
assignopt( {'time limits', 'tlim'}, 'tlim', v );
assignbopt( 'byrate', 'byrate', v );
assignbopt( '3d', 'use_3d', v );

if doblack, colordef black; end

if byrate
  [foo perm] = sort( cellfun(@length, st), 'descend' );
end

if isempty(perm)
  perm = 1:length(st);
end

nclr = size(clr,1);
st = st(perm);

% restrict to time limits
if ~isempty(tlim)
  f = @(x,y)(x(x >= y(1) & x <= y(2)));
  [y{1:length(st)}] = deal(tlim);
  st = cellfun(f, st, y, 'uniformoutput', 0);
end
  
for ii = 1:length(st)
  clrid = mod(ii - 1, nclr) + 1;
  nt = length(st{ii});
  u = repmat(st{ii}(:)', 3, 1);
  v = repmat([1-th/2; 1+th/2; nan], 1, nt) + spc*(ii-1);
  if use_3d
      plot3(u(:), v(:), v(:)*0,  '-', 'color', clr(clrid,:), 'linewidth', linewidth); 
      view(2);
      hold on
  else
    plot(u(:), v(:), '-', 'color', clr(clrid,:), 'linewidth', linewidth);
    hold on
  end
end
if ~isempty(tlim), xlim(tlim); end

if ~isempty( ii ), ylim(spc*[0 ii]+.5); end
hold off

if doblack, colordef white; end

set(pan,'Motion','horizontal','Enable','on')
set(zoom,'Motion','horizontal','Enable','on')
set(gca,'xtick',[],'ytick',[]);
set(gca,'ydir','Normal')
