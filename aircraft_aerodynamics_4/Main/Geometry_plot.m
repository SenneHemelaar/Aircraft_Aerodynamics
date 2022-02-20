clear; close all;

lw    = 1.68;
cwr   = 1.6;
lamw  = 0.4;
cwt   = lamw*cwr;
Lamw  = deg2rad(30);
phiw_list = linspace(0, 90, 7);

%%% Main Wing %%%
X_wing = [0 3.5 5.5 5.5];
Y_wing = [0 14 14 0];
Z_wing = [0 0 0 0];

%%% Winglet %%%
X_winglet = [X_wing(3),...
             X_wing(3) - cwr,...
             3.5 + (lw * cos(Lamw)),...
             3.5 + (lw * cos(Lamw)) + cwt];

for i = 1:7
phiw  = deg2rad(phiw_list(i));

Y_winglet{i} = [14,...
             14,...
             14 + (lw * sin(phiw)),...
             14 + (lw * sin(phiw))];

Z_winglet{i} = [0,...
             0,...
             lw*cos(phiw),...
             lw*cos(phiw)];
end

figure(1)
hold on;
view(3)
patch(X_wing,Y_wing,Z_wing,'black')
patch(X_wing,Y_wing,Z_wing,'black')
patch(X_winglet,Y_winglet{1},Z_winglet{1},'red')
patch(X_winglet,Y_winglet{2},Z_winglet{2},'blue')
patch(X_winglet,Y_winglet{3},Z_winglet{3},'green')
patch(X_winglet,Y_winglet{4},Z_winglet{4},'yellow')
patch(X_winglet,Y_winglet{5},Z_winglet{5},'red')
patch(X_winglet,Y_winglet{6},Z_winglet{6},'blue')
patch(X_winglet,Y_winglet{7},Z_winglet{7},'green')
axis('equal')
ylim([12 16]);


