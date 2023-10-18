close all
clear all
clc

%number of blades
N = 2;

%diameter of the propeller
dia=1.6;

%tip radius
r_tip = dia/2;

%start blade
r_hub = 0.1*R;     %percentage of tip radius

%number of sections
rstep=(r_tip-r_hub)/10; %number of steps 
r1 = r_hub:rstep:r_tip;

%interpolated formula of chord length
chord = (20.672*r1^6 - 84.504*r1^5 + 124.14*r1^4 - 82.258*r1^3 + 21.475*r1^2 + 0.1965*r1 + 0.3403)*2;

%pitch distance in meters.
pitch=1.0;

%engine speed in RPM
RPM=2100;
%thickness to chord ratio for propeller section (constant with radius)
tonc=0.12*chord;
%standard sea level atmosphere density
rho=1.225;
%RPM --> revs per sec
n=RPM/60.0;
%rps --> rads per sec
omega=n*2.0*pi;
% use 10 blade segments (starting at 10% R (hub) to R)
xs=0.1*R;
xt=R;
rstep=(xt-xs)/10;
r1=xs:rstep:xt;


%calculate results for a range of velocities from 1 to 60m/s
for V_inf=1:60
    %initialise sums
    thrust=0.0;
    torque=0.0;
    %loop over each blade element
    for j=1:size(r1,2)
        r=r1(j);
        %calculate local blade element setting angle
        theta=atan(pitch/2/pi/r); %beta
        %guess initial values of inflow and swirl factor
        a_i_0=0.1;
        b_i_0=0.01;
        %set logical variable to control iteration
        finished=false;
        %set iteration count and check flag
        sum=1;
        itercheck=0;
        fail = 0;
        while (~finished)
            %%%====Equations 3 & 4====%%%
            
            %%% Calculate velocities
            %axial velocity
            Vax=V_inf*(1+a_i_0);
            %disk plane velocity
            Vtan=omega*r*(1-b_i_0);
            
            % Equation 3
            %local velocity at blade
            Vp=sqrt(Vax^2+Vtan^2);
            
            % Equation 4
            %flow angle %which flow angle
            phi=atan2(Vax,Vtan);
            %blade angle of attack 
            alpha=theta-phi;
            
            %%%====Equations 1 & 2 ====%%%
            %%% Get lift and drag coefficients
            % lift coefficient
            cl=6.2*alpha; %Too simplified lift coefficent
            %drag coefficient
            cd=0.008-0.003*cl+0.01*cl*cl;  %Too simplified drag coefficient
            
            
            % Equation 1 divided by dr
            DTdr = 0.5 * rho * chord * Vp^2 * (cl*cos(phi) - cd*sin(phi)) * N;
            
            % Equation 2 divided by dr
            DQdr = 0.5 * rho * chord * Vp^2 * (cd*cos(phi) + cl*sin(phi)) * N*r;

            %%%====Equations 5 & 6====%%%
            %inflow and swirl
            a_i_1 = DTdr/(rho * 4 * pi * r * V_inf^2 * (1+a_i_0));
            b_i_1 = DQdr/(rho * 4 * pi * r^3 * V_inf * (1+a_i_0) * omega);
            
            
            %stabilise iteration
            a_i_1 = 0.5 * a_i_0 + 0.5 * a_i_1; 
            b_i_1 = 0.5 * b_i_0 + 0.5 * b_i_1;
            
            %check for convergence
            if (abs(a_i_1-a_i_0)<1.0e-5) && (abs(b_i_1-b_i_0)<1.0e-5)
                finished=true;
            end


            a_i_0=a_i_1;
            b_i_0=b_i_1;
            %increment iteration count
            sum=sum+1;
            %check to see if iteration stuck
            if (sum>500)
                finished=true;
                itercheck=1;
                fail = fail + 1;
            end

        end
    thrust=thrust+DTdr*rstep;
    torque=torque+DQdr*rstep;
        
    end
    CT(V_inf)=thrust/(rho*n^2*dia^4);
    CP(V_inf)=torque/(rho*n^2*dia^5);
    J(V_inf)=V_inf/(n*dia);
    alpha_lib(V_inf) = torque;
    if CT(V_inf) < 0
        eff(V_inf) = 0;
    else
        eff(V_inf)=J(V_inf)/2.0/pi*(CT(V_inf)/CP(V_inf));
    end
end

Jmax=max(J);
CTmax=max(CT);

f1 = figure(1);
plot(J,CT,J,CP);
title('Thrust and Torque Coefficients')
xlabel('Advance Ratio (J)');
ylabel('Ct, Cq');
legend('Ct','Cq');
axis([0 Jmax 0 1.1*CTmax ]);

f2 = figure(2);
plot(J,eff);
title('Propeller Efficiency');
xlabel('Advance Ratio (J)');
ylabel('Efficiency');
axis([0 Jmax 0 1 ]);