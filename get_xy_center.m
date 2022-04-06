function dat_JB=get_xy_center(dat_JB)
%get the x and y for each timepoint for each larva
for i=1:length(dat_JB.xspine)
    dat_JB.x{i,1}=dat_JB.xspine{i,1}(:,6);
    dat_JB.y{i,1}=dat_JB.yspine{i,1}(:,6);
end 
end 