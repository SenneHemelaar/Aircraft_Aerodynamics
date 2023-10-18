function [thrust, power, torque, Vax_list,...
        Vtan_list, r_R_list, thrust_list, ...
        torque_list, a_list, b_list] = ...
        SolveForFreeStreamVelocity(V_inf, sec, dr, ...
        tipR, rootR, pitch, omega, N, n, r_steps, rho, ...
        Cl_a, Cl_a_sec_fn, Cd_a, Cd_a_sec_fn, chord, ...
        PrandtlCorrection, DiagnosticInfo)
    
    TSR = n*tipR/V_inf;
    %initialise sums
    thrust=0.0;
    torque=0.0;
    power =0.0;
    %loop over each blade element
    Vax_list  = zeros(length(r_steps),1);
    Vtan_list = zeros(length(r_steps),1);
    a_list    = zeros(length(r_steps),1);
    b_list    = zeros(length(r_steps),1);
    for j=1:r_steps
        
        r = (sec(j) + 0.5*dr(j))*tipR;
        r_R = r/tipR;
        r_R_list(j) = r_R;
        theta = interp1(pitch(1,:),pitch(2,:),sec(j),'spline');
        theta = atan(theta/2/pi/r);
        a_i_0=0.1;
        b_i_0=0.01;
        %set logical variable to control iteration
        finished=false;
        %set iteration count and check flag
        sum_iter=1;
        itercheck=0;
        
        while (~finished)
            %axial velocity
            Vax=V_inf*(1+a_i_0);
            Vax_list(j) = Vax;
            
            %disk plane velocity
            Vtan=omega*r*(1-b_i_0);
            Vtan_list(j) = Vtan;
            
            %local velocity at blade
            Vp=sqrt(Vax^2+Vtan^2);
            
            %flow angle %which flow angle
            phi=atan2(Vax,Vtan);
            
            %blade angle of attack 
            alpha=theta-phi;
            
            %%%====Equations 1 & 2 ====%%%
            %%% Get lift and drag coefficients
            % lift coefficient
%             if j == 1
%                 cl = 0.01;
%             elseif j == 2
%                 cl = 0.045*alpha;
%             else
                cl = interp1(Cl_a.(Cl_a_sec_fn{j})(1,:), ... 
                             Cl_a.(Cl_a_sec_fn{j})(2,:), ...
                             rad2deg(alpha));
                if isnan(cl)
                    cl = interp1(Cl_a.(Cl_a_sec_fn{j})(1,:), ... 
                             Cl_a.(Cl_a_sec_fn{j})(2,:), ...
                             rad2deg(alpha), 'linear','extrap');
                end
%             end
            
            % drag coefficient
            if j == 1
                cd = 0.4;
            elseif j == 2
                cd = 0.1;
            else 
                cd = 1/(interp1(Cd_a.(Cd_a_sec_fn{j})(1,:), ... 
                                Cd_a.(Cd_a_sec_fn{j})(2,:), ...
                                rad2deg(alpha),'spline')/cl);
            end
            
            % chord
            chord_j = interp1(chord(1,:), ... 
                              chord(2,:), ...
                              sec(j),'spline');
                          
            % Equation 1 divided by dr
            DTdr = 0.5 * rho * chord_j * Vp^2 * (cl*cos(phi) - cd*sin(phi)) * N;
            
            % Equation 2 divided by dr
            DQdr = 0.5 * rho * chord_j * Vp^2 * (cd*cos(phi) + cl*sin(phi)) * N*r;

            %%%====Equations 5 & 6====%%%
            %inflow and swirl
            a_i_1 = DTdr/(rho * 4 * pi * r * V_inf^2 * (1+a_i_0));
            b_i_1 = DQdr/(rho * 4 * pi * r^3 * V_inf * (1+a_i_0) * omega);
            
            %%% Prandtl tip and root losses %%% 
            for xx=1
            if PrandtlCorrection
                temp1 = -N/2 * ((tipR/tipR)-r_R)/r_R *sqrt( 1+ ((TSR*r_R)^2)/((1-a_i_1)^2));
                F_tip = 2/pi * acos(exp(temp1));
                F_tip(isnan(F_tip)) = 0;

%                 temp2 = N/2 * ((rootR/tipR)-r_R)/r_R *sqrt( 1+ ((TSR*r_R)^2)/((1-a_i_1)^2));
%                 F_root = 2/pi * acos(exp(temp2));            
%                 F_root(isnan(F_root)) = 0;

                pr_corr =  F_tip;%*F_root;
            else
                pr_corr = 1;
            end
            
            if pr_corr < 1e-4
                pr_corr = 1e-4;
            end
            end
            
            %stabilise iteration
            a_i_1 = 0.5 * a_i_0 + 0.5 * a_i_1*pr_corr; 
            b_i_1 = 0.5 * b_i_0 + 0.5 * b_i_1*pr_corr;
            
            %check for convergence
            if (abs(a_i_1-a_i_0)<1.0e-5) && (abs(b_i_1-b_i_0)<1.0e-5)
                finished=true;
            end

            a_i_0=a_i_1;
            b_i_0=b_i_1;
            %increment iteration count
            sum_iter=sum_iter+1;
            %check to see if iteration stuck
            if (sum_iter>50)
                finished=true;
                itercheck=1;
                if DiagnosticInfo
                    fprintf('NaN for V_inf = %.2f \n',V_inf)
                end
            end
            
        end

        thrust_list(j) = DTdr*dr(j)*tipR;
        torque_list(j) = DQdr*dr(j)*tipR;
        a_list(j)      = a_i_0;
        b_list(j)      = b_i_0;
        
        thrust= thrust + DTdr*dr(j)*tipR;
        torque= torque + DQdr*dr(j)*tipR;
        power = power +  DQdr*dr(j)*tipR*omega;
        
        if DiagnosticInfo
            fprintf('j = %i |  cl = %.2f, cd = %.2f, r =  %.2f, dr =  %.2f,',[j cl cd r dr(j)*tipR]);
            fprintf('chord = %.2f, pitch = %.2f, alpha = %.2f, Ftotal = %f\n',[chord_j,rad2deg(theta),rad2deg(alpha),F_tip*F_root]);
        end
    end
end
    