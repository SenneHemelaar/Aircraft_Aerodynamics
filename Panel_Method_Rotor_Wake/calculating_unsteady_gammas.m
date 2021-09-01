function [gamma_b_i, gamma_w] = calculating_unsteady_gammas(v_g,N,u_0,w_0, t_n, gamma_shedded, x_v_w, z_v_w)
% this function calculates the vortex strength on each thin airfoil panel
% and the stength of the shedded vortex. 

% defining empty lists
u = zeros(N+1);             % x component induced vel. by bound vortex j at control point i 
w = zeros(N+1);             % z component induced vol. by bound vortex j at control point i
a = zeros(N+1);             % amtric coefficients: dot product induced vel. with normal vector at cp
RHS = zeros(1,N+1);         % right hand side of linear equation 
u_w = zeros(1, N);          % x compontnt induced vel. by wake vortex k at control point i 
w_w = zeros(1, N);          % z component induced vel. by wake vortex k at control point i


% adding new shedded vortex coordinate to vortex coordinate list for loop
x_v = [v_g.x_v, v_g.x_v(end) + 0.25 * (x_v_w(1) - v_g.x_v(end))];
z_v = [v_g.z_v, v_g.z_v(end) + 0.25 * (z_v_w(1) - v_g.z_v(end))];

% The wake is continous and is discretized as the center of it's range
wake_center_x = x_v_w - 0.5 * (x_v_w - [v_g.x_v(end), x_v_w(1:end-1)]);
wake_center_z = z_v_w - 0.5 * (z_v_w - [v_g.z_v(end), z_v_w(1:end-1)]);

for i=1:N+1 % begin loop over control points + total circulation conserved
    for j=1:N+1 % begin loop over bound vortices + first vortices 
        
        % total circulation should be conserved by assigning last row 
        % values 1. corrseponding RHS term is assigned total previously 
        % calculated vortex trength
        if i == N+1 
            a(i,j) = 1;
        else
            [u(i,j),w(i,j)]= VOR2D(v_g.x_c(i),v_g.z_c(i),x_v(j),z_v(j), 1); % unit strength 1 is assigned to Circulation at this point 
            a(i,j)=u(i,j)*v_g.x_n(i)+w(i,j)*v_g.z_n(i);
        end
    end % end loop over bound vortex + wake vortex
    
    % calculating induced velocities by wake vorices
    if i == N+1
        RHS(i) = -sum(gamma_shedded);
    else
        for k = 1:t_n-1 % loop over wake vorices
            [u_ind, w_ind] = VOR2D(v_g.x_c(i),v_g.z_c(i), wake_center_x(k+1), wake_center_z(k+1), gamma_shedded(k));
            u_w(i) = u_w(i) + u_ind;
            w_w(i) = w_w(i) + w_ind; 
        end 
        RHS(i)=-(u_0*v_g.x_n(i)+w_0*v_g.z_n(i))-(u_w(i)*v_g.x_n(i)+w_w(i)*v_g.z_n(i));
    end 
    

end % end loop over contol points

RHSt = transpose(RHS); %The vector needs to be transposed from column to rows to solve the linear equation.

gammas = linsolve(a,RHSt); 


gamma_b_i = gammas(1:end-1);
gamma_w = gammas(end);

end