function [X_wing, Y_wing, Z_wing, X_winglet, Y_winglet, Z_winglet]...
         = Geometry_plot_function(lw,  Lam,  phiw,  cwr,  lam,  epsR,  epsT,...
                                  cR,  cT,   L,     NC1,  NC2,  NS1,   NS2)

%%%======================Preliminary Calculations=======================%%%
cwR     = cwr*cT;       % Winglet root chord length
cwT     = lam*cwR;      % Winglet tip chord length [m]
Lw      = 2*L*lw;       % Total length of winglet [m]

%%%=============Main Wing and Winglet Coordinates Generation============%%%
%   Main wing geometry
wing_geom    = [ 0.0 ,0.0,0.0,cR,0.0 ;
                cR-cT, L ,0.0,cT,0.0];

%   Winglet geometry             
winglet_geom = [       cR-cwR      ,       L      ,     0.0    ,cwR,-epsR ;
                cR-cwR+Lw*tand(Lam),L+Lw*sind(phiw),Lw*cosd(phiw),cwT,-epsT];

%   Zero Matrices for Plotting the Mesh
X_wing    = zeros(NC1,NS1);
Y_wing    = zeros(NC1,NS1);
Z_wing    = zeros(NC1,NS1);

LE_X_wing = cosspace(wing_geom(1,1),wing_geom(2,1),NS1);
LE_Y_wing = cosspace(wing_geom(1,2),wing_geom(2,2),NS1);

X_winglet = zeros(NC2,NS2);
Y_winglet = zeros(NC2,NS2);
Z_winglet = zeros(NC2,NS2);

if winglet_geom(2,1) > winglet_geom(1,1)
    LE_X_winglet = cosspace(winglet_geom(1,1),...
        winglet_geom(2,1),NS2);
elseif winglet_geom(1,1) > winglet_geom(2,1)
    LE_X_winglet = fliplr(cosspace(winglet_geom(2,1),...
        winglet_geom(1,1),NS2));
else
    LE_X_winglet = ones(1,NS2)*winglet_geom(1,1);
end

%   Define all points at the winglet trailing edge at each section in X

if winglet_geom(2,1)+winglet_geom(2,4)>winglet_geom(1,1)+winglet_geom(1,4)
    TE_X_winglet = cosspace(winglet_geom(1,1)+winglet_geom(1,4),...
        winglet_geom(2,1)+winglet_geom(2,4),NS2);
elseif winglet_geom(2,1)+winglet_geom(2,4)<...
        winglet_geom(1,1)+winglet_geom(1,4)
    TE_X_winglet = fliplr(cosspace(winglet_geom(2,1)+...
        winglet_geom(2,4),winglet_geom(1,1)+winglet_geom(1,4),NS2));
else
    TE_X_winglet = ones(1,NS2)*(winglet_geom(2,1)+winglet_geom(2,4));
end

%   Define all points at the winglet leading edge at each section in Y

if winglet_geom(2,2) > winglet_geom(1,2)
    LE_Y_winglet = cosspace(winglet_geom(1,2),winglet_geom(2,2),NS2);
elseif winglet_geom(2,2) < winglet_geom(1,2)
    LE_Y_winglet = fliplr(cosspace(winglet_geom(2,2),...
        winglet_geom(1,2),NS2));
else
    LE_Y_winglet = ones(1,NS2)*(winglet_geom(1,2));
end

%   Define all points at the winglet leading edge at each section in Z
if winglet_geom(2,3) > winglet_geom(1,3)
    LE_Z_winglet = cosspace(winglet_geom(1,3),winglet_geom(2,3),NS2);
elseif winglet_geom(2,3) < winglet_geom(1,3)
    LE_Z_winglet = fliplr(cosspace(winglet_geom(2,3),...
        winglet_geom(1,3),NS2));
else
    LE_Z_winglet = ones(1,NS2)*winglet_geom(1,3);
end

%   Finish defining all coordinates
for i = 1:NS1
    tempX = cosspace(LE_X_wing(i),wing_geom(1,4),NC1);
    tempY = LE_Y_wing(i);
    for j = 1:NC1
        X_wing(j,i) = tempX(j);
        Y_wing(j,i) = tempY;
    end
end

for i = 1:NS2
    tempX = cosspace(LE_X_winglet(i),TE_X_winglet(i),NC2);
    tempY = LE_Y_winglet(i);
    tempZ = LE_Z_winglet(i);
    for j = 1:NC2
        X_winglet(j,i) = tempX(j);
        Y_winglet(j,i) = tempY;
        Z_winglet(j,i) = tempZ;
    end
end

end

