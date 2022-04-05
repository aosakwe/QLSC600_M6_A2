function deg=calculate_heading_direction(xspine,yspine,vector)
%this function just calculates the heading direction at each timepoint for
%each larva based on JB data
%vector=[1 0] or [-1 0] for the different position; [-1 0] for odor in the
%left side 
%in the code you can put

deg={};
for i=1:length(xspine)
    x1=xspine{i,1}(:,1);
    x2=xspine{i,1}(:,7);
    x3=x1-x2;
    
    y1=yspine{i,1}(:,1);
    y2=yspine{i,1}(:,7);
    y3=y1-y2;
    % change (-1,0) based on the position of odor
    degree=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));
    deg{i,1}=0-degree;
    clear degree x1 x2 x3 y1 y2 y3
end
end