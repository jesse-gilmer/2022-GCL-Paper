function outlist = assignopt( key, lhsvarname, inlist )
% key - name of the option to retrieve
% lhsvarname - name of variable to assign value to, in caller frame
% inlist - option list
%
% Example:
% assignopt( {'buffer size', 'buff'}, 'x', v )
% will assign the value of 'buff' to x in the calling frame, iff
% there is a 'buff' option specified in the list v.

outlist = inlist;
if( getopt( key, inlist ) )
    [val outlist] = getopt( key, inlist, 'value' );
    assignin( 'caller', lhsvarname, val );
end



