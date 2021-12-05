function helpFun(name)
%
% Help for functions in Command Window
%
% Syntax
%     helpFun('myFunction')
%
% Description
%     Display the help text for the function specified by NAME. This
%     function can be used as a workaround if the native Matlab help
%     function does not show properly the help section in the command
%     window, and displays all successively commented lines immediately
%     following the function declaration.
%
% Input arguments
%     NAME [char(1 x :inf)] is the name of the function for which help is
%     requested.
%
% Output arguments
% None
%
% Example
%     % Show the help section of the function NNO.m
%     helpFun('NNO')
%
% Copyright (c) 2021 by George Papazafeiropoulos
% _________________________________________________________________________
%


name = strtrim(name);
% Add file extension if needed
if isempty(regexp(name,'.m$','once'))
    name = [name,'.m'];
end
fid = fopen(name);
% Search for the line containing the 'function' keyword
foundfunc = false;
while ~foundfunc
    C = textscan(fid,'%s',1,'delimiter','\n');
    L = C{1}{1};
    foundfunc = ~isempty(regexp(strtrim(L),'^function','once'));
end
% Switch showing if the first line beginning with '%' has been found 
sw1=false;
% Scan all lines after the function definition
while true
    C = textscan(fid,'%s',1,'delimiter','\n');
    L = C{1}{1};
    percentpos = regexp(L,'^%');
    if ~isempty(percentpos)
        % If line contains '%', then print it to command window
        fprintf('%s\n',L(percentpos+1:end));
        sw1=true;
    elseif sw1
        % If line does not contain '%' and at least one line beginning with
        % '%' has already been found, then stop printing
        break
    end
end
fclose(fid);
end

