% Program to solve non-Newtonian polymer flow
% 1D non-Newtonian flow
% Matthew T. Balhoff
% February 18, 2017

close all; clear all; clc

% Read in the input file
[reservoir fluid numerical well P BC]= inputfile;

%Create a vicosity plot by looping through shear rates 
figure (1)
shear_rate = logspace(-3,4, 100);
for i= 1:length(shear_rate)
    viscosity (i) = fluid.visc_inf +(fluid.visc_p-fluid.visc_inf)...
        *(1+(fluid.lamda*shear_rate(i))^fluid.a)^((fluid.n-1)/fluid.a);   
end
loglog(shear_rate,viscosity)
ylabel ('apparent viscosity (cp)')
xlabel ('apparent shear rate (1/s)')
print -djpeg -r300 'visc_curve.jpg'

% Implicitly calculate reservoir pressure through time
reservoir.vel= zeros(numerical.N,1);reservoir.shear= zeros(numerical.N,1);
time=0; n = 1;                  % Initialize the time as zero;  % "n" is the time step (an integer) 
P_plot=zeros(numerical.N,ceil(numerical.t_final/numerical.dt));
BHP=zeros(200,2);
while time < numerical.t_final            % Continue while time is less than final time (defined above)
   P_plot(:,n)=P; visc_plot(:,n)=fluid.visc; vel_plot(:,n)=reservoir.vel; shear_plot(:,n)=reservoir.shear;
   P_old=P;                     % Create a "dummy" vector as placeholder for previous Pressure 
    
   if numerical.iterate== 'explicit'        % explicit updates to the noinlinear problem
       % Compute block velcoity, viscosity, and shear rate
        [fluid, reservoir] = visc_iterate (fluid, reservoir, numerical,well,P);
        
       % Cpompute the arrays with the new viscosity
        [T,B,Q,jprod] =myarrays(reservoir, fluid, numerical, well, P, BC);
        
        % Solve  the preessure field
        P = (T+B/numerical.dt)\(B*P_old/numerical.dt + Q);
   elseif numerical.iterate== 'picard  '
        [fluid, reservoir] = visc_iterate (fluid, reservoir, numerical,well,P);
        [T,B,Q,jprod] =myarrays(reservoir, fluid, numerical, well, P, BC);
        P = (T+B/numerical.dt)\(B*P_old/numerical.dt + Q);
        error=1.0; tol=1E-6;
        while error > tol
            Ptemp=P;
            [fluid, reservoir] = visc_iterate(fluid, reservoir, numerical,well,P);
            [T,B,Q,jprod] =myarrays(reservoir, fluid, numerical, well, P, BC);
            P = (T+B/numerical.dt)\(B*P_old/numerical.dt + Q);
            error= sqrt(abs(P-Ptemp)'*abs(P-Ptemp));
        end
   elseif numerical.iterate == 'newton  '
        error=1.0; tol=1E-6;
        [J, F] =myJacobian(reservoir, fluid, numerical, well, P, BC, P_old);
        while error > tol
            P = P - J\F;
            [fluid, reservoir] = visc_iterate (fluid, reservoir, numerical,well,P);
            [J, F] =myJacobian(reservoir, fluid, numerical, well, P, BC, P_old);
        
            error= sqrt(F'*F);
            [T,B,Q, jprod] =myarrays(reservoir, fluid, numerical, well, P, BC);
        end
   end
   n=n+1;                           % Update time level by one time step.
   time = time + numerical.dt;      % update actual time by adding the time increment, dt
end

% make the plots
postprocess(P_plot,visc_plot, vel_plot, shear_plot,numerical,reservoir,fluid);

%time=1:n-1;
%plotyy(time, BHP(:,1), time, BHP(:,2))