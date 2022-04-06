function []=plot_specific_trajectory(list,x,y,turn_x,turn_y,t0_idx,t1_idx)
%list is the index of trajectory that you want to plot, it will just be the row number
figure(1)
set(gcf,'Name','Trajectory');
hold on
for i=1:length(list)

    plot(x{list(i),1},y{list(i),1},'g.','MarkerSize',10);
    hold on
    plot(x{list(i),1}(1),y{list(i),1}(1),'k.','MarkerSize',18);

    if ~isempty(turn_x)&& ~isempty(turn_y)
        idx=[];
        for j=1:length(t0_idx{list(i),1})
            idx=[idx;(t0_idx{list(i),1}(j):t1_idx{list(i),1}(j))'];
        end
        plot(x{list(i),1}(idx),y{list(i),1}(idx),'c.','MarkerSize',10);
        plot(turn_x{list(i),1},turn_y{list(i),1},'r.','MarkerSize',15);
        plot(turn_x{list(i),2},turn_y{list(i),2},'b.','MarkerSize',15);

    end

end
legend({'run','trajectory start','turn','turn start','turn stop'});
end