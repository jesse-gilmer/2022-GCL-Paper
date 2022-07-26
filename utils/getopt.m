function [opt outlist] = getopt(key, inlist, varargin)

%GETOPT   Get option from list.
%
%   TF = GETOPT(KEY, INLIST) returns 0 or 1 depending on whether
%   the keyword KEY is present in the cell INLIST.
%
%   [TF OUTLIST] = GETOPT(...) also returns a copy of the INLIST 
%   but with the keyword KEY removed if it were present.
%
%   VAL = GETOPT( ... , 'value') instead returns the element of
%   INLIST immediately following the keyword or an empty array
%   if the keyword was not found.
%
%   [VAL OUTLIST] = GETOPT( ... , 'value') also removes the keyword
%   value pair from the input list and returns it in OUTLIST.

%   Author: E.V.Lubenov, Caltech
%   Date: 31-Oct-2005

% copy the inlist

% check args
error(nargchk(2,3,nargin));

% check if returning boolean or value
if (nargin == 3) && strcmpi(varargin{1}, 'value')
  vflg = 1;
else
  vflg = 0;
end

% copy inlinst
outlist = inlist;

% cell of keys is used to handle synonyms
if ~iscell(key), key = {key}; end

% ensure that keyword occurs only once in inlist
bool = zeros(size(inlist));
for ii = 1:length(key), bool = bool | strcmpi(inlist, key{ii}); end
if sum(bool) > 1, 
  error([mfilename ': more than one keyword "' key{1} '" encountered.']);
end

kind = find(bool);
if vflg
  if ~isempty(kind)
    % make sure key is followed by a values
    if (kind + 1 > length(inlist))
      error([mfilename ': no value following keyword "' key{1} '".']);
    end
    opt = inlist{kind+1};
    outlist(kind:kind+1) = [];
  else
    opt = [];
  end
else
  opt = ~isempty(kind);
  outlist(kind) = [];
end
