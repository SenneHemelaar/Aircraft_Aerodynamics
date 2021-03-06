function [gamma] = calculating_gamma(v_g, N, U_inf, W_inf, n_x, n_z)

x_c = v_g.x_c;
z_c = v_g.z_c;
x_v = v_g.x_v;
z_v = v_g.z_v;


u=zeros(1,N);
w=zeros(1,N);
a=zeros(1,N);
RHS=zeros(1,N);

for i=1:N % begin loop over control points  
    for j=1:N % begin loop over bound vortices  
        [u(i,j),w(i,j)]= VOR2D(x_c(i),z_c(i),x_v(j),z_v(j), 1);
       
        
        a(i,j)=u(i,j)*n_x(i)+w(i,j)*n_z(i);
 
    end
    
    RHS(i)=-(U_inf*n_x(i)+W_inf*n_z(i));
  
end

RHSt=transpose(RHS); %The vector needs to be transposed from column to rows to solve the linear equation.

gamma=linsolve(a,RHSt); 

end

