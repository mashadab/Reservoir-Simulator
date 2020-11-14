function [T,k_half,velocity_half,shear_half, visc_half] = Thalf(i,j, reservoir, fluid, numerical,P)

if P(i)>P(j)
    visc =fluid.visc(i);
else
    visc = fluid.visc(j);
end

T=(reservoir.perm(i)*reservoir.h*reservoir.W)/(numerical.dx(i)*visc);


% k_half= (numerical.dx(i)+numerical.dx(j))/(numerical.dx(i)/reservoir.perm(i) + numerical.dx(j)/reservoir.perm(j));
% dx_half = (numerical.dx(i)+numerical.dx(j))/2;
% 
% % iterate to determine interblock viscosity
% visc_half = 0.5*fluid.visc(i) + 0.5*fluid.visc(j);
% visc_half=fluid.visc_p;
% phi_half=(reservoir.phi(i)+reservoir.phi(j))/2;
% error=1.0;
% tol=1E-6;
% while error> tol
%     visc_old=visc_half;
%     velocity_half = (k_half/visc_half)*((P(i)-P(j))/dx_half)*6.33E-03;
%     shear_half = 112.3*6.0*((3*fluid.n +1)/(4*fluid.n))^(fluid.n/(fluid.n-1))...
%         *4*velocity_half/(sqrt(8*k_half*phi_half));
%     %shear_half = 112.3*6.0*((3*fluid.n +1)/(4*fluid.n))^(fluid.n/(fluid.n-1))...
%     %    *velocity_half/(sqrt(k_half*phi_half)); 
%     visc_half = fluid.visc_inf +(fluid.visc_p-fluid.visc_inf)...
%         *(1+(fluid.lamda*shear_half)^fluid.a)^((fluid.n-1)/fluid.a);
%     error= abs((visc_half-visc_old)/visc_half);
% end
% T=k_half*(reservoir.h*reservoir.W)/(visc_half*dx_half);
% 
% visc_half;