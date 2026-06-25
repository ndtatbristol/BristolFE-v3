clear
close all

%Switch to current folder so that files can be found and created in correct
%place
cd(fileparts(mfilename('fullpath')));


readme_fname = '..\readme.md';
function_summary_heading = '## Function summary';
fn_summary_line = 'SUMMARY';

lines = readlines(readme_fname);
idx = find(contains(lines, function_summary_heading), 1, 'first');
if ~isempty(idx)
    writelines(lines(1:idx-1), readme_fname);
end

writelines([function_summary_heading, ""], readme_fname, WriteMode = "append");

function_fnames = dir('fn*.m');
for i = 1:numel(function_fnames)
    [~, nm, ~] = fileparts(function_fnames(i).name);
    
    [usage, author, summary, inputs, outputs] = fn_parse_fn_file(function_fnames(i).name);

    if summary ~= ""
        lines = ["### " + nm,
            "",
            summary,
            "```",
            usage,
            "```"];
        writelines(lines, readme_fname, WriteMode = "append");
    else
        fprintf([nm, ' missing SUMMARY\n']);
    end
end

