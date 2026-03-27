function fn_console_output(short_comment, varargin)
global LOG_TO_FILE
if numel(varargin) < 1 || isempty(varargin{1})
    long_comment = short_comment;
else
    long_comment = varargin{1};
end
if numel(varargin) < 2
    indent = 1;
else
    indent = varargin{2};
end
comment_indent_level = fn_get_comment_indent_level;
if indent
    indent_str = repmat('  - ', 1, comment_indent_level);
else
    indent_str = '';
end
% fprintf(['Indent level %i, ', short_comment],  comment_indent_level);
switch fn_get_comment_verbosity
    case 'low'
        str = [indent_str, short_comment];
    case 'high'
        str = [indent_str, long_comment];
    otherwise
end
fprintf(str);
if exist('LOG_TO_FILE', "var")
    try
        fid = fopen(LOG_TO_FILE,"a");
        fprintf(fid, str);
        fclose(fid);
    catch
        % warning(['Could not write to file ', LOG_TO_FILE])
    end
end
end