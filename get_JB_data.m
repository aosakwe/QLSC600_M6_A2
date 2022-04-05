function [dat,uname]=get_JB_data(folder_name,ii)
% get the data from JB analysis
%create empty strcture where we store the data
choredir =folder_name;
dat.xspine = {};
dat.yspine ={};
dat.AN={};
dat.et={};

%get all the subfolder name 
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
%cellfun: apply a function to the array to get the last section of the name of the file 
end

if ispc == 1
    splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
else
    splits = cellfun(@(x) [x{end-1},'/',x{end}], splits, 'UniformOutput', false);
end
[uname,na,nb] = unique(splits);%find how many conditions we have based on the filename 

idx = find(nb == ii);
for jj = idx'%(1:3)'% for demostraton we just want to load three folders, if you want to load all, delete '(1:3)'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*trx.mat']));%open the trx file in each folder
    f = contains({fname.name},['trx.']);
    fname = fullfile(fname(f).folder,fname(f).name);
    fileID = open(fname);
    dat.et= vertcat(dat.et, {fileID.trx_t.t}');
    dat.xspine = vertcat (dat.xspine, {fileID.trx_t.x_spine}');%save the specific data into each field of the structure 
    dat.yspine = vertcat (dat.yspine, {fileID.trx_t.y_spine}');
    dat.AN= vertcat(dat.AN, {fileID.trx_t.numero_larva_num}');
    clear fileID

end
fn=fieldnames(dat);

for i=1:length(fn)
    dat.(fn{i})=vertcat(dat.(fn{i}){:});
end
end