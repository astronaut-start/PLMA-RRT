function [result, scene, options] = run_planner(algorithm, environment, seed, overrides)
if nargin<1, algorithm='PLMA'; end
if nargin<2, environment='S1'; end
if nargin<3, seed=1; end
if nargin<4, overrides=struct; end
id=upper(string(algorithm));
ids=["GB","BI","F","Q","BAI","IBIP","BIAPF","PLMA"];
assert(isscalar(id)&&any(id==ids),'plma:Algorithm', ...
    'Use GB, BI, F, Q, BAI, IBIP, BIAPF, or PLMA.');
validateattributes(seed,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',2^32-1});
assert(isstruct(overrides)&&isscalar(overrides),'plma:Options','overrides must be a scalar struct.');
root=fileparts(mfilename('fullpath'));
source=fullfile(root,'src'); native=fullfile(source,'plma_repaired');
addpath(source,native,fullfile(native,'common'),fullfile(native,'controller'), ...
    fullfile(root,'environments'),'-begin');
functions={'planners.gb_karaman2011','planners.bi_jordan2013','planners.frrt_liao2021', ...
    'planners.qrrt_jeong2019','planners.bai_xu2025','planners.ibip_fan2023', ...
    'planners.biapf_huang2024','planners.plma_repaired'};
for name=[functions,{'rr.collision','rr.obstacle_hit','plma_adaptive','collision'}]
    actual=which(name{1});
    assert(startsWith(lower(actual),[lower(source) filesep]), ...
        'plma:PathConflict','%s resolves outside this repository: %s',name{1},actual);
end
scene=load_environment(environment);
options=planner_options(char(id),scene);
fields=fieldnames(overrides);
if id=="PLMA"
    assert(all(ismember(fields,{'max_iter','time_limit'})), ...
        'plma:Options','The PLMA adapter accepts only max_iter and time_limit overrides.');
end
for k=1:numel(fields), options.(fields{k})=overrides.(fields{k}); end
validateattributes(options.max_iter,{'numeric'},{'scalar','real','finite','integer','positive'});
validateattributes(options.time_limit,{'numeric'},{'scalar','real','positive','nonnan'});
previous=rng; restore=onCleanup(@()rng(previous)); %#ok<NASGU>
rng(double(seed),'twister');
planner=str2func(functions{find(id==ids,1)});
timer=tic;
result=planner(scene,options);
result.wall_time=toc(timer);
result.algorithm_id=char(id);
result.environment_id=sprintf('S%d',scene.scene);
result.seed=double(seed);
end
