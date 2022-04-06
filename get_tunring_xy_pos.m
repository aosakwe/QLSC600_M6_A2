function [turn_x,turn_y]=get_tunring_xy_pos(x,y,t0_idx,t1_idx)

for i=1:length(x)
    idx0=t0_idx{i,1}; idx1=t1_idx{i,1};
    if isempty(idx0)||isempty(idx1)
        turn_x{i,1}=[];turn_y{i,1}=[];
    else
        x_t0=x{i,1}(idx0);% x pos for the startpoint of each turning event
        x_t1=x{i,1}(idx1);% x pos for the endpoint of each turning event
        x_mid=x_t0+((x_t1-x_t0)./2);% x pos for the midpoint of each turning event
        turn_x{i,1}=x_t0;
        turn_x{i,2}=x_t1;
        turn_x{i,3}=x_mid;

        y_t0=y{i,1}(idx0);% x pos for the startpoint of each turning event
        y_t1=y{i,1}(idx1);% x pos for the endpoint of each turning event
        y_mid=y_t0+((y_t1-y_t0)./2);% x pos for the midpoint of each turning event
        turn_y{i,1}=y_t0;
        turn_y{i,2}=y_t1;
        turn_y{i,3}=y_mid;

        end
end

end