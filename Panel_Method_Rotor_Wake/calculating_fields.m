function [abs_vel_field, u_vel_field,w_vel_field,pres_field] = calculating_fields(gamma, x_v, z_v, x_field, z_field, rho, u_0, w_0)
% this function finds the pressure and absolute velocities from the 

% velocities induced by the bound vortices on the airfoil panels 
% The velocities are calculated on a pre defined grid (x-field)x(z-field)


abs_vel_field = zeros(length(x_field), length(z_field));
u_vel_field = zeros(length(x_field), length(z_field)) + u_0;
w_vel_field = zeros(length(x_field), length(z_field)) + w_0;
pres_field = zeros(length(x_field), length(z_field));

for i = 1:length(z_field)
    for j = 1:length(x_field)
        for k =1:length(gamma(:,1))
            [u,w] = VOR2D(x_field(j),z_field(i),x_v(k),z_v(k), gamma(k));

            
            u_vel_field(i,j) = u_vel_field(i,j) + u;
            w_vel_field(i,j) = w_vel_field(i,j) + w;
        end

        abs_vel_field(i,j) = (u_vel_field(i,j)^2+w_vel_field(i,j)^2)^0.5;
        pres_field(i,j) = (0.5*rho*(u_0^2+w_0^2) -0.5*rho*abs_vel_field(i,j)^2)/(0.5*rho*(u_0^2+w_0^2)) ;
    end 

end 
end 

    