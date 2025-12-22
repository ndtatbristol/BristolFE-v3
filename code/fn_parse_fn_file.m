function [usage, author, summary, inputs, outputs] = fn_parse_fn_file(fname)
%USAGE
%   [usage, summary, inputs, outputs] = fn_parse_fn_file(fname)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Parses function file and extracts contents of each part of help
%INPUTS
%   fname - file name of function to parse
%OUTPUTS
%   usage, author, summary, inputs, outputs - as per names. Each one
%   returned as a single string with whitespace etc. stripped out.
%--------------------------------------------------------------------------

flines = readlines(fname);

usage = fn_extract(flines, 'USAGE');
author = fn_extract(flines, 'AUTHOR');
summary = fn_extract(flines, 'SUMMARY');
inputs = fn_extract(flines, 'INPUTS');
outputs = fn_extract(flines, 'OUTPUTS');

end

function content = fn_extract(flines, heading)
idx = find(contains(flines, heading), 1, 'first');
if ~isempty(idx)
    content = join(flines(idx+1:end));
    content = strtrim(regexprep(erase(extractBetween(content, 1, regexp(content, '%(?=\S)', 'once')), "%"), '\s+', ' '));
else
    content = "";
end

end

