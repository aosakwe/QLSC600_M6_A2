function [t0_idx,t1_idx,t0s,t1s]=find_ind_t0s_t1s(et,t0s,t1s)
% this function will find the corresponding index for the t0s and t1s
%et is time from JB, t0s and t1s is from JAABA
for i=1:length(et)
    dat0=t0s{i,1};
    dat1=t1s{i,1};
    time=et{i,1};
     %find the t0s and t1s that is out of the time range and delete them
    ind0=find(~ismember(dat0,time))';
    ind1=find(~ismember(dat1,time))';
    ind=[ind0,ind1];
    dat0(ind)=[];
    dat1(ind)=[];

    t0s{i,1}(ind)=[];
    t1s{i,1}(ind)=[];
    %get the index for t0s and t1s based on et
    idx0=find(ismember(time,dat0));
    idx1=find(ismember(time,dat1));
    t0_idx{i,1}=idx0;
    t1_idx{i,1}=idx1;
end

clearvars -except  t*
end 