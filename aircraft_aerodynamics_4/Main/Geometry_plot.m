%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                Geometry Plot for AVL in Matlab                      %%%
%%%                     Author: Senne Hemelaar                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%%%==========================INPUTS GUIDE===============================%%%
% lw      : Winglet Length                          [% span]
% Lam     : Winglet Leading Edge Sweep              [degrees]
% phiw    : Winglet Cant Angle                      [degrees]
% cwr     : Winglet Root Chord                      [% main wing chord]
% lam     : Winglet Taper Ratio                     [cwt/cwr]
% epsR    : Winglet Root Twist Angle                [degrees]
% epsT    : Winglet Tip Twist Angle                 [degrees]
% cR      : Main Wing Root Chord                    [m]
% cT      : Main Wing Tip Chord                     [m]
% L       : Main Wing Half-Span                     [m]
% NC1     : No. of Chordwise Sections Main Wing     [-]
% NC2     : No. of Chordwise Sections Winglet       [-]
% NS1     : No. of Spanwise Sections Main Wing      [-]
% NS2     : No. of Spanwise Sections Winglet        [-]

%%%=============================INPUTS==================================%%%
lw        = 0.6;
Lam       = 8.2046;
phiw      = 0;
cwr       = 0.9072;
lam       = 0.5307;
epsR      = 0;
epsT      = 0;
cR        = 5.5;
cT        = 2.0;
L         = 14;
NC1       = 10;
NC2       = 8;
NS1       = 22;
NS2       = 22;

color_string = {'#2B3A67'; '#5C6D70'; '#A37774'; '#E88873'; '#E0AC9D';...
                '#2B3A67'; '#5C6D70'; '#A37774'; '#E88873'; '#E0AC9D';...
                '#2B3A67'; '#5C6D70'; '#A37774'; '#E88873'; '#E0AC9D'};
% '#484A47' : Black/Grey   ----   '#C490D1' : Purple

%%%=====================CALL GEOMETRY FUNCTION==========================%%%
[X_wing, Y_wing, Z_wing, X_winglet, Y_winglet, Z_winglet]....
= Geometry_plot_function(lw,  Lam,  phiw,  cwr,  lam,  epsR,  epsT,...
                         cR,   cT,     L,  NC1,  NC2,  NS1,   NS2 );

%%%============================PLOTTING=================================%%%
figure(1)
view(300,30)
hold on; grid on;
axis equal
mesh(X_wing,Y_wing,Z_wing,'FaceAlpha',1.0,'EdgeColor','black',...
    'FaceColor',[1 0 1],'LineWidth',0.75);
% mesh(X_wing,-Y_wing,Z_wing,'FaceAlpha',1.0,'EdgeColor','black',...
%     'FaceColor','#C490D1','LineWidth',0.75);
mesh(X_winglet,Y_winglet,Z_winglet,'FaceAlpha',1.0,'EdgeColor',...
    'black','FaceColor',[1 0 1],'LineWidth',0.75);
% mesh(X_winglet,-Y_winglet,Z_winglet,'FaceAlpha',1.0,'EdgeColor',...
%     'black','FaceColor','#C490D1','LineWidth',0.75);
xlabel('$X$ [m]','interpreter','latex');
% ylim([10 16]); xlim([2 6.5]);zlim([0 0.2])
ylabel('$Y$ [m]','interpreter','latex');
zlabel('$Z$ [m]','interpreter','latex');
set(gcf,'position',[100,100,1500,750])















% lw    = 1.68;
% cwr   = 1.6;
% lamw  = 0.4;
% cwt   = lamw*cwr;
% Lamw  = deg2rad(30);
% phiw_list = linspace(0, 90, 7);

% %%% Main Wing %%%
% X_wing = [0 3.5 5.5 5.5];
% Y_wing = [0 14 14 0];
% Z_wing = [0 0 0 0];
% 
% %%% Winglet %%%
% X_winglet = [X_wing(3),...
%              X_wing(3) - cwr,...
%              3.5 + (lw * cos(Lamw)),...
%              3.5 + (lw * cos(Lamw)) + cwt];
% 
% for i = 1:7
% phiw  = deg2rad(phiw_list(i));
% 
% Y_winglet{i} = [14,...
%              14,...
%              14 + (lw * sin(phiw)),...
%              14 + (lw * sin(phiw))];
% 
% Z_winglet{i} = [0,...
%              0,...
%              lw*cos(phiw),...
%              lw*cos(phiw)];
% end
% 
% figure(1)
% hold on; grid on;
% view(3)
% patch(X_wing,Y_wing,Z_wing,'black')
% patch(X_wing,-Y_wing,Z_wing,'black')
% patch(X_winglet,Y_winglet{1},Z_winglet{1},'red')
% patch(X_winglet,-Y_winglet{1},Z_winglet{1},'red')
% % patch(X_winglet,Y_winglet{2},Z_winglet{2},'blue')
% % patch(X_winglet,Y_winglet{3},Z_winglet{3},'green')
% % patch(X_winglet,Y_winglet{4},Z_winglet{4},'yellow')
% % patch(X_winglet,Y_winglet{5},Z_winglet{5},'red')
% % patch(X_winglet,Y_winglet{6},Z_winglet{6},'blue')
% % patch(X_winglet,Y_winglet{7},Z_winglet{7},'green')
% xlabel('$x$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% zlabel('$z$','interpreter','latex')
% axis('equal')
% ylim([-16 16]);
% zlim([0 2]);

