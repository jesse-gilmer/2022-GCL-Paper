function outlist = assignbopt( key, lhsvarname, inlist )
% function outlist = assignbopt( key, lhsvarname, inlist )
% Assign binary option (0=false, 1=true)
% key - name of the option to retrieve
% lhsvarname - name of variable to assign value to, in caller frame
% inlist - option list
%
% Example:
% assignopt( {'verbose', 'v'}, 'verb', varargin )
% will assign the value of 1 to 'verb' in the calling frame, iff
% there is a 'verbose' option specified in the list varargin.
%
% Returns the input list, with the first instance of the keyword
% removed, if found.

outlist = inlist;
if( getopt( key, inlist ) )
  assignin( 'caller', lhsvarname, 1 );
  loc = find(strcmp(key,outlist));
  outlist = outlist([ 1:(loc(1)-1) (loc(1)+1):length(outlist)]);
end



