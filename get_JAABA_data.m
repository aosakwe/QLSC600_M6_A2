function [dat,uname]=get_JAABA_data(folder_name,ii)
%create empty strcture with 3 fields 
choredir=folder_name;
dat.t0s={};
dat.t1s={};
dat.AN={};
%get all the subfolders inside this directory 
full = dir(choredir);
n = 1;
sear = fullfile(choredir,'*');
d = 1;
while length(d) > 0
    sear = fullfile(sear,'*');
    d = dir(sear);
    full = vertcat(full,d);
    n = n+1;
end

filt = [full.isdir];
full = full(filt);
names = {full.name}';

expr = '^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
filt = regexp(names,expr);
filt = cellfun(@(x) ~isempty(x), filt);
d = full(filt);

names = {d.folder}';

if ispc == 1
    splits = cellfun(@(x) split(x,'\'), names, 'UniformOutput', false);
else
    splits = cellfun(@(x) split(x,'/'), names, 'UniformOutput', false);
end %cellfun: apply a function to the array
if ispc == 1
    splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
else
    splits = cellfun(@(x) [x{end-1},'/',x{end}], splits, 'UniformOutput', false);
end

[uname,na,nb] = unique(splits);

idx = find(nb == ii);

for jj = idx(1:3)'% for demostraton we just want to load three folders, if you want to load all, delete '(1:3)'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*scores_turn.mat']));%load the table and save it to the structure 
    f = contains({fname.name},['scores_turn.']);
    fname = fullfile(fname(f).folder,fname(f).name);

    fileID = open(fname);
    dat.t0s = vertcat (dat.t0s, {fileID.score.turn_starts_in_seconds}');
    dat.t1s = vertcat (dat.t1s, {fileID.score.turn_stops_in_seconds}');
    dat.AN =vertcat (dat.AN, {fileID.score.animal_num}');
    clear fileID fname f
end
clearvars -except uname dat
end