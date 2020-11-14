function [fluid reservoir] = visc_iterate (fluid, reservoir, numerical,well,P)


% loop through all other blocks
for i=1:numerical.N-1
   fluid.visc_old(i) = fluid.visc(i);
   
   tol=1E-6; error =1;
   while error>tol
        reservoir.vel(i) = abs(6.33E-03*Thalf (i,i+1,reservoir, fluid, numerical,P)*abs(P(i)-P(i+1)))/(reservoir.h*reservoir.W);
        
        reservoir.shear(i) = 112.3*6.0*((3*fluid.n +1)/(4*fluid.n))^(fluid.n/(fluid.n-1))...
        *4*reservoir.vel(i)/(sqrt(8*reservoir.perm(i)*reservoir.phi(i)));
    
        fluid.visc(i) = fluid.visc_inf +(fluid.visc_p-fluid.visc_inf)...
        *(1+(fluid.lamda*reservoir.shear(i))^fluid.a)^((fluid.n-1)/fluid.a);
        
        error = abs(fluid.visc(i)-fluid.visc_old(i));
        fluid.visc_old(i) = fluid.visc(i);
   end
end

% Last block is just Q/A;
reservoir.vel(numerical.N) = abs(well.rates(2)/(reservoir.h*reservoir.W));
reservoir.shear(numerical.N) = 112.3*6.0*((3*fluid.n +1)/(4*fluid.n))^(fluid.n/(fluid.n-1))...
        *4*reservoir.vel(numerical.N)/(sqrt(8*reservoir.perm(i)*reservoir.phi(i)));
fluid.visc(numerical.N) = fluid.visc_inf +(fluid.visc_p-fluid.visc_inf)...
        *(1+(fluid.lamda*reservoir.shear(numerical.N))^fluid.a)^((fluid.n-1)/fluid.a);
