function [T,B,Q,jprod] =myarrays(reservoir, fluid, numerical, well, P, BC)

% Set up matrix T and vector Q
T = sparse(numerical.N,numerical.N); B=sparse(numerical.N,numerical.N); 
Q = sparse(numerical.N,1);
for i=1:numerical.N
   
    if i==1
        T(i,i+1) = -Thalf(i,i+1, reservoir, fluid,numerical,P); 
        if strcmp(BC.type(1,:),'neumann  ') 
            T(i,i) = T(i,i)-T(i,i+1);
        elseif strcmp(BC.type(1,:),'dirichlet') 
            T(i,i)=  T(i,i) - T(i,i+1) + 2*Thalf(i,i,reservoir, fluid,numerical,P); 
            Q(i,1) = 2*Thalf(i,i,reservoir, fluid,numerical,P)*BC.value(1)*6.33E-03;
        end
        
    elseif i==numerical.N
        T(i,i-1) = -Thalf(i,i-1, reservoir, fluid,numerical,P);  
         if strcmp(BC.type(2,:),'neumann  ')           
            T(i,i) = T(i,i)-T(i,i-1);
         elseif strcmp(BC.type(2,:),'dirichlet') 
            T(i,i)=  T(i,i) - T(i,i-1) + 2*Thalf(i,i,reservoir, fluid,numerical,P); 
            Q(i,1) = 2*Thalf(i,i,reservoir, fluid,numerical,P)*BC.value(2)*6.33E-03;
        end
    else        
        T(i,i-1)= -Thalf(i,i-1,reservoir, fluid,numerical,P);            % Upper diagonal
        T(i,i+1)= -Thalf(i,i+1,reservoir, fluid,numerical,P);            % Lower diagonal
        T(i,i) = - (T(i,i-1) + T(i,i+1));   % Main diagonal
    end
    B(i,i) = (reservoir.h*reservoir.W)*numerical.dx(i)*reservoir.phi(i)*fluid.cf;
end

for j=1:length(well.grids)
    i=well.grids(j);
    Q(i)=Q(i) + well.Q(i);
    jprod(j) = 6.33E-03*(2*pi*reservoir.perm(i)*reservoir.h)/(fluid.visc_inf.*log(0.2.*numerical.dx(i))./well.radii(j));
end

%G= -0.433*6.33E-03*T*z;
T=6.33E-03*T;