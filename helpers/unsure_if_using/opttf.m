function tf = opttf(pval,okarg,mfile)
%OPTTF determines whether input options are true or false

% Copyright 2003-2007 The MathWorks, Inc.
% $Revision: 1.3.4.3 $   $Date: 2007/09/11 11:42:46 $


if islogical(pval)
    tf = all(pval);
    return
end
if isnumeric(pval)
    tf = all(pval~=0);
    return
end
if ischar(pval)
    truevals = {'true','yes','on','t'};
    k = any(strcmpi(pval,truevals));
    if k
        tf = true;
        return
    end
    falsevals = {'false','no','off','f'};
    k = any(strcmpi(pval,falsevals));
    if k
        tf = false;
        return
    end
end
if nargin == 1
    % return empty if unknown value
    tf = logical([]);
else
    okarg(1) = upper(okarg(1));
    xcptn = MException(sprintf('Bioinfo:%s:%sOptionNotLogical',mfile,okarg),...
        '%s must be a logical value, true or false.',...
        upper(okarg));
    xcptn.throwAsCaller;
end