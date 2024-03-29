function [a , a_prime, r_R_annulus, f_normal, f_tangential, gamma, alpha, inflowangle] = solveStreamtube(U_0, rho, r1_R, r2_R, rootradius_R, tipradius_R , omega, R, n_blades, chord, twist, alfaMat, clMat, cdMat, yaw,  max_iterations, error_limit, iteratingDesignParameters)
    %% design optimalization parameter
    iterationsize = 1.01;
    
    %% dimensional parameters
    r_R_annulus = (r2_R + r1_R)/2;  %middle-point of annulus
    Area = pi*((r2_R*R)^2-(r1_R*R)^2);
    dr = r2_R - r1_R;
    
    %% tracking variables
    chord_initial = chord;
    twist_initial = twist;
    last_edited = 'chord';
    previous_C_p = 0;

    a = 0; 
    a_prime = 0;
    a_initial = 0;

    %% Fixing a for convergence
    if iteratingDesignParameters
         Ct_fixed = 0.75;
         a_goal = induction_factor(Ct_fixed);
    end
    
    %% Actual iteration process    
    for j = 1:max_iterations
        if(yaw == 0 || j == 1)
            U_rotor = U_0 * (1-a);
            U_tangential = (1+a_prime)*omega*r_R_annulus*R;
            [f_normal, f_tangential, gamma, alpha, inflowangle] = loadBladeElement(U_rotor, U_tangential, r_R_annulus, chord, twist, alfaMat, clMat, cdMat);
        else
            [f_normal, f_tangential, gamma, alpha, inflowangle] = solveAzimuthalElement(U_0, R, r_R_annulus, chord, twist, alfaMat, clMat, cdMat, a, a_prime, omega, yaw);
        end
        loadaxial = f_normal*R*dr*n_blades;
        CT_initial = loadaxial/(0.5*Area*U_0^2);
        CT = CT_initial;
        CP = dr * f_tangential * r_R_annulus * R * n_blades * omega/(0.5*U_0^3*pi*R^2*rho);
        anew = induction_factor(CT);
        
        %{
        %}

        Prandtl = PrandtlTipRootCorrection(r_R_annulus, rootradius_R, tipradius_R, omega*R/U_0, n_blades, anew);
        %to ensure we are not dividing by zero
        if (Prandtl < 0.0001)
            Prandtl = 0.0001; 
        end
        
        
        anew = anew/Prandtl;
        a = 0.75*a+0.25*anew;

        a_prime = f_tangential*n_blades/(2*pi*U_0*(1-a)*omega*2*(r_R_annulus*R)^2);
        a_prime = a_prime/Prandtl;
        
        if (j == 1)
            a_initial = a;
        end
        
        % test convergence of solution, by checking convergence of axial induction
        if(j == max_iterations)
            disp('max iterations reached');
            disp('didnt actually converge');
        end
        
        if (abs(a-anew) < error_limit && (~iteratingDesignParameters || a > 0.5) || iteratingDesignParameters && abs(a_goal-anew) < error_limit)
            disp("Did converge")
            break
        end
        
        if (iteratingDesignParameters)
            [chord, twist, last_edited] = design_optimalization(last_edited, CP, previous_C_p, chord, twist, iterationsize, twist_initial);
            previous_C_p = CP;
        end
    end
