%Author: Lucas Philipp
%run code_for_JB_JAABA_analysis.m first

%odor is on the right x=0
%origin is top right

% Dat_JAABA – detection of turning events
% t0s : beginning time of turn
% t1s : end time of turn
% AN: animal ID
% t0_idx : Index (flame number) of beginning of turn
% t1_idx : index (flame number) of end of turn
% pre_deg: trajectory angle of before turning
% post_deg: trajectory angle of after turning
% turn_x: x position of turn cols: (start, end, midpoint)
% turn_y: y position of turn cols: (start, end, midpoint)
 
% Dat_JB
% {xspine} : x coordinate of spine 
% {yspine} : y coordinate of spine
% AN: animal ID
% Et: time
% X : X coordinate of center of mass
% Y: X coordinate of center of mass

%plot head direction vs turning frequency
%turning frequency could be defined in two ways
%1. time spent turning with a certain pre_deg/total time (not just turn time) for all trajectories
%2. total # of turns with a certain pre_deg/total # of turns

pol_bins=100; %number of bins in angle space (same for theta and phi)
sqbins = 25; %number of square bins (must divide into 250)

%2. is easier
fns = fieldnames(dat_JAABA);
t0s=dat_JAABA.(fns{1});
t1s=dat_JAABA.(fns{2});
AN_JAABA=dat_JAABA.(fns{3});
t0_idx=dat_JAABA.(fns{4});
t1_idx=dat_JAABA.(fns{5});
pre_deg=dat_JAABA.(fns{6});
post_deg=dat_JAABA.(fns{7});
x_turn=dat_JAABA.(fns{8});
y_turn=dat_JAABA.(fns{9});

fns2 = fieldnames(dat_JB);
xspine=dat_JB.(fns2{1});
yspine=dat_JB.(fns2{2});
AN_JB=dat_JB.(fns2{3});
Et=dat_JB.(fns2{4});
X=dat_JB.(fns2{5});
Y=dat_JB.(fns2{6}); 

pre_deg_all=[];
post_deg_all=[];
t0s_all=[];
t1s_all=[];
x_turn_all=[];
y_turn_all=[];
X_all=[];
Y_all=[];

tot_time = 0;
for i=1:size(pre_deg,1)
    pre_deg_all = [pre_deg_all; pre_deg{i}];
    post_deg_all = [post_deg_all; post_deg{i}];
    t0s_all = [t0s_all; t0s{i}];
    t1s_all = [t1s_all; t1s{i}];
    x_turn_all = [x_turn_all; x_turn{i}]; %this keeps only the x positions of the start  of the turn
    y_turn_all = [y_turn_all; y_turn{i}]; %this keeps only the y positions of the start  of the turn
    tot_time = tot_time + Et{i}(end);
    X_all = [X_all; X{i}];
    Y_all = [Y_all; Y{i}];
end

figure
polarhistogram(pre_deg_all,pol_bins)
title({'head direction v.s. turning frequency=total # of turns with a certain start direction','(all trajectories)/total # of turns (all trajectories)'})
%1.
dur=zeros(size(pre_deg,1),1);
dur=t1s_all-t0s_all;
%get indexes of turns in each angle bin
pol_edges=linspace(-180,180,pol_bins+1);
time_pre_deg=zeros(1,pol_bins);
for i =1:pol_bins
    ind=find(pre_deg_all>=pol_edges(i) & pre_deg_all<=pol_edges(i+1));
    time_pre_deg(i)=sum(dur(ind));
end

levels=zeros(1,pol_bins);
for i = 1:pol_bins
    levels(i)=(pol_edges(i)+pol_edges(i+1))/2;
end
figure
polarplot([levels*pi/180 levels(1)*pi/180], [time_pre_deg./tot_time time_pre_deg(1)./tot_time])
title({'head direction v.s. turning frequency==time spent turning with a certain start direction','(all trajectories)/total time spent turning (all trajectories)'})

%position vs turning frequency
%time spent turning in this grid square/total time spent in this grid square

%divide 250x250 into grid squares of size 10
xedges = linspace(0,250,sqbins+1);
yedges = linspace(0,250,sqbins+1);

%take position as the position at the start of the turn
%contatenate x_turn_all and y_turn_all for find()
pos_turn_all=zeros(size(x_turn_all,1),2);
pos_turn_all(:,1)=x_turn_all;
pos_turn_all(:,2)=y_turn_all;

%time spend turning in the grid square
time_turn_xy = zeros(250/sqbins);

for i = 1:sqbins
    for j = 1:sqbins
    ind=find(pos_turn_all(:,1)>=xedges(i) & pos_turn_all(:,1)<=xedges(i+1) & pos_turn_all(:,2)>=yedges(j) & pos_turn_all(:,2)<=yedges(j+1));
    time_turn_xy(i,j)=sum(dur(ind));
    end
end

figure('color','w');
imagesc(time_turn_xy); % This will plot a 2D matrix with value of each weight encoded by a colour.
ax=gca;
colorbar; % The color bar will tell you what value each colour corresponds to.
set(ax,'DataAspectRatio',[1,1,1]); % Specific every axis to have the same ratio
% Alternatively, you can use the command "axis square"
res=20; %resoliution of the colour map
CMnegative = [linspace(0,1,res)',linspace(0,1,res)',ones(res,1)]; % a res by 3 array ranging from [0,0,1] to [1,1,1]
CMpositive = [ones(res,1),linspace(1,0,res)',linspace(1,0,res)']; % a res by 3 array ranging from [1,1,1] to [1,0,0]
CM = [CMnegative;CMpositive(2:end,:)];
colormap(CM);
ax.CLim = [0 max(max(time_turn_xy))];

set(gca,'XTick',linspace(0,25,sqbins))
set(gca,'YTick',linspace(0,25,sqbins))

xtickcell={};
ytickcell={};
sto=linspace(0,250,sqbins+1);
stoinv=linspace(250,0,sqbins+1);
for i=1:sqbins
xtickcell(end+1) = {num2str(stoinv(i))}; %origin is top right
ytickcell(end+1) = {num2str(sto(i))};
end
set(gca,'yticklabel',ytickcell)
set(gca,'xticklabel',xtickcell)

title('Total time [s] spent turning in each grid cell');
xlabel('$x$ position','interpreter','latex')
ylabel('$y$ position','interpreter','latex')

%time spend in the grid square
time_grid = zeros(250/sqbins);

time_spent=[];
for i=1:size(pre_deg,1) 
time_spent = [time_spent; [0; diff(Et{i})]];
end

XY_all=zeros(size(X_all,1),2);
%contatenate arrays for find()
XY_all(:,1)=X_all;
XY_all(:,2)=Y_all;

for i = 1:sqbins
    for j = 1:sqbins
    ind=find(XY_all(:,1)>=xedges(i) & XY_all(:,1)<=xedges(i+1) & XY_all(:,2)>=yedges(j) & XY_all(:,2)<=yedges(j+1));
    time_grid(i,j)=sum(time_spent(ind));
    end
end

figure('color','w');
imagesc(time_grid); % This will plot a 2D matrix with value of each weight encoded by a colour.
ax=gca;
colorbar; % The color bar will tell you what value each colour corresponds to.
set(ax,'DataAspectRatio',[1,1,1]); % Specific every axis to have the same ratio
% Alternatively, you can use the command "axis square"
colormap(CM);
ax.CLim = [0 max(max(time_grid))];

set(gca,'XTick',linspace(0,25,sqbins))
set(gca,'YTick',linspace(0,25,sqbins))
set(gca,'yticklabel',ytickcell)
set(gca,'xticklabel',xtickcell)

title('Total time [s] spent in each grid cell');
xlabel('$x$ position','interpreter','latex')
ylabel('$y$ position','interpreter','latex')

figure('color','w');
imagesc(time_turn_xy./time_grid); % This will plot a 2D matrix with value of each weight encoded by a colour.
ax=gca;
colorbar; % The color bar will tell you what value each colour corresponds to.
set(ax,'DataAspectRatio',[1,1,1]); % Specific every axis to have the same ratio
% Alternatively, you can use the command "axis square"
colormap(CM);
ax.CLim = [0 max(max(time_turn_xy./time_grid))];

set(gca,'XTick',linspace(0,25,sqbins))
set(gca,'YTick',linspace(0,25,sqbins))
set(gca,'yticklabel',ytickcell)
set(gca,'xticklabel',xtickcell)

title('time spent turning in this grid square/total time spent in this grid square');
xlabel('$x$ position','interpreter','latex')
ylabel('$y$ position','interpreter','latex')

div=time_turn_xy./time_grid;
[row, col] = find(div >=1);
diag(time_grid(row,col))
diag(time_turn_xy(row,col))

%head direction vs turning direction
%see: https://www.mathworks.com/matlabcentral/answers/6625-3-d-histogram-in-spherical-coordinates

%contatenate for find()
turn_dir_all=zeros(size(pre_deg_all,1),2);    
turn_dir_all(:,1)=pre_deg_all;
turn_dir_all(:,2)=post_deg_all;

turn_pair = zeros(pol_bins);

for i = 1:pol_bins
    for j = 1:pol_bins
    ind=find(turn_dir_all(:,1)>=pol_edges(i) & turn_dir_all(:,1)<=pol_edges(i+1) & turn_dir_all(:,2)>=pol_edges(j) & turn_dir_all(:,2)<=pol_edges(j+1));
    turn_pair(i,j)=size(ind,1);
    end
end

%H is your histogram data
H = turn_pair;

%%Set up the figure
colordef(figure,'black');
theta_vec = linspace(0,2*pi,pol_bins);
phi_vec = linspace(0,2*pi,pol_bins);
[theta,phi] = meshgrid(theta_vec,phi_vec);
Hmax = max(H(:));
r = 0.03*Hmax; %Box size
polar(nan,max(max(H.*cos(phi))));
hold all;
%%Make the Histogram
for kk = 1:numel(theta_vec);
    for jj = 1:numel(phi_vec);
        X=r*([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5);
        Y=r*([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5);
        Z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*H(jj,kk);
        h= patch(X,Y,Z,0*X+H(jj,kk),'edgecolor','none');
          rotate(h,[0 0 1],45,[0 0 0]);
          rotate(h,[0 1 0],90 - 180/pi*phi_vec(jj),[0 0 0]);
          rotate(h,[0 0 1],180/pi*theta_vec(kk),[0 0 0]);
      end;
  end;
%%Adjust the plot
[Xs,Ys,Zs] = sphere(size(theta,2)+1);
hs = surf(Hmax*Xs,Hmax*Ys,Hmax*Zs);
set(hs,'facecolor','none','edgecolor','w','edgealpha',0.2)
camlight;
set(gca,{'xtick' 'ytick' 'ztick' 'vis' 'clim'},{[] [] [] 'on' [0 Hmax]});
axis equal vis3d;
box on;
view(3);
colorbar
drawnow;
