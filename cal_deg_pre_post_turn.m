function dat_JAABA=cal_deg_pre_post_turn(vector,dat_JB,dat_JAABA)

%the pre and post turning heading direction
%vector is the odor position which will be [-1 0], you can check this by
%plotting the x and y position 
for i=1:length(dat_JAABA.t0_idx)
    idx0=dat_JAABA.t0_idx{i,1};
    idx1=dat_JAABA.t1_idx{i,1};
    if isempty(idx0)||isempty(idx1)
        %if there is no turning event detect in here, go to the next animal
        dat_JAABA.post_deg{i,1}=[];dat_JAABA.pre_deg{i,1}=[];
        continue
    end
    % calculate heading directio using x and y 5 points before the turning event
    ind=find(idx0==1|idx1==length(dat_JB.x{i,1}));%if a turning event happens in the first and the last point of the trajectory, delete it
    idx0(ind)=[];idx1(ind)=[];
    dat_JAABA.t0s{i,1}(ind)=[]; dat_JAABA.t1s{i,1}(ind)=[]; dat_JAABA.t0_idx{i,1}(ind)=[]; dat_JAABA.t1_idx{i,1}(ind)=[];
    if isempty(idx0)||isempty(idx1)
        %if there is no turning event detect here, delete the turning event
        %detect in JAABA as this turning event happens in the start and end point
        %because it is hard to calculate the pre-turning direction
        dat_JAABA.t0s{i,1}=[]; dat_JAABA.t1s{i,1}=[]; dat_JAABA.t0_idx{i,1}=[]; dat_JAABA.t1_idx{i,1}=[];
        dat_JAABA.post_deg{i,1}=[];dat_JAABA.pre_deg{i,1}=[];
        continue
    end

    %% calculate the pre_deg
    idx0_pre=idx0-5;
    neg_value=find(idx0_pre<1); %if there is a negative value of index, change it to 1, so that we can calculate the pre-deg
    idx0_pre(neg_value)=1;

    x1=dat_JB.x{i,1}(idx0);% then use the x and y value to get a vector 
    y1=dat_JB.y{i,1}(idx0);

    x2=dat_JB.x{i,1}(idx0_pre);
    y2=dat_JB.y{i,1}(idx0_pre);

    x3=x1-x2;
    y3=y1-y2;
    
    pre_deg=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));% calculate the heading direction
    clear x* y*
    dat_JAABA.pre_deg{i,1}=0-pre_deg;
    %% calculate the post_deg 
    idx1_post=idx1+5;
    ind=find(idx1_post>length(dat_JB.x{i,1})); %if the index is longer than the total length, change it 
    idx1_post(ind)=length(dat_JB.x{i,1});

    x2=dat_JB.x{i,1}(idx1);% then use the x and y value to get a vector 
    y2=dat_JB.y{i,1}(idx1);

    x1=dat_JB.x{i,1}(idx1_post);
    y1=dat_JB.y{i,1}(idx1_post);

    x3=x1-x2;
    y3=y1-y2;
    post_deg=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));% calculate the heading direction
    dat_JAABA.post_deg{i,1}=0-post_deg;
end

end
