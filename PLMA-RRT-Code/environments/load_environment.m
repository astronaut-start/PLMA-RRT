function s = load_environment(environment)
if isnumeric(environment)
    validateattributes(environment,{'numeric'},{'scalar','integer','>=',1,'<=',3});
    id=double(environment);
else
    value=upper(string(environment));
    assert(isscalar(value)&&any(value==["S1","S2","S3"]), ...
        'plma:Environment','Use S1, S2, S3, or numeric 1, 2, 3.');
    id=str2double(extractAfter(value,1));
end
file=fullfile(fileparts(mfilename('fullpath')),sprintf('S%d.json',id));
s=jsondecode(fileread(file));
s.start=double(s.start(:)'); s.goal=double(s.goal(:)');
if isnumeric(s.obstacles), s.obstacles=num2cell(s.obstacles,2); end
s.obstacles=cellfun(@(v)double(v(:)'),s.obstacles,'UniformOutput',false);
assert(isequal(size(s.space),[3 2])&&numel(s.start)==3&&numel(s.goal)==3, ...
    'plma:Environment','Expected a three-dimensional scene.');
end
