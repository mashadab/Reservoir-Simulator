function [reservoir fluid numerical well P BC]= inputfile

% Define reservoir properties 
reservoir.L = 4000;            % Reservoir length (ft)
reservoir.h = 20;               % reservoir thickness (h)
reservoir.W = 500;             % area (ft2)
reservoir.perm = 200*ones(1,50);            % permeability (mD)
reservoir.phi = 0.25*ones(1,50);            % porosity

% Define fluid properties
fluid.cf = 1E-5;                % fluid compressibility (1/psi)
fluid.visc_inf = 1.0;           % high shear viscosity (cp)
fluid.visc_p = 100.0;            % high shear viscosity (cp)
fluid.n = 0.5;                  % shear thinning index (0<n<1)
fluid.lamda = 50;
fluid.a= 2.0;                
fluid.visc = fluid.visc_p*ones(length(reservoir.perm),1);

% Define numerical parameters for discretized solution
numerical.N=length(reservoir.perm); % Number of grids          
numerical.dx = (reservoir.L/numerical.N)*ones(numerical.N,1);             % Grid size (ft)
numerical.dt=0.1;               % Timestep (days)
numerical.t_final = 10;        % final time (days)      
numerical.method = 'implicit';  % implicit, explicit, or mixed_CN
numerical.iterate = 'picard  '; %explicit, picard , or Newton   

% compute positions of block centers
numerical.x(1)=numerical.dx(1)/2;
for i=2:numerical.N
    numerical.x(i)=numerical.x(i-1) + 0.5*(numerical.dx(i-1)+numerical.dx(i));
end

% Well data
well.locations = [1 3999];          % position, feet
well.rates = [1000  -1000];       % rate, ft3/day or BHP, psi
well.radii = [0.25 0.25];           % well radius, ft
well.type = [1 1];                  % 1= rate well; 2= BHP well
for i=1:length(well.locations)
    for j=1:numerical.N
        if well.locations(i) < numerical.x(j)+numerical.dx(j)/2 && well.locations(i) > numerical.x(j)-numerical.dx(j)/2
            well.grids(i)=j;
        end
    end
end
well.Q(well.grids)=well.rates;      % vector of rates, ft3/day

BC.type = ['neumann  ';'neumann  '];
BC.value =[0.0; 0.0]; % value of boundary condition, ft3/day or psi
P= 3000*ones(numerical.N,1);