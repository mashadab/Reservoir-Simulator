function [J,F] =myJacobian(reservoir, fluid, numerical, well, P, BC, P_old)

Q = sparse(numerical.N,1);
for j=1:length(well.grids)
    i=well.grids(j);
    Q(i)=Q(i) + well.Q(i);
    jprod(j) = 6.33E-03*(2*pi*reservoir.perm(i)*reservoir.h)/(fluid.visc_inf.*log(0.2.*numerical.dx(i))./well.radii(j));
end

epsilon=1.0E-6;
% Set up matrix J and vector F
J = sparse(numerical.N,numerical.N); B=sparse(numerical.N,numerical.N); 
F = sparse(numerical.N,1); T = sparse(numerical.N,numerical.N);
for i=1:numerical.N
    B(i,i) = (reservoir.h*reservoir.W)*numerical.dx(i)*reservoir.phi(i)*fluid.cf;
    if i==1
        T(i,i+1) = -6.33E-03*Thalf(i,i+1, reservoir, fluid,numerical,P); 
        F(i,1) = T(i,i+1)*(P(i+1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        
        Ptemp=P;
        P(i+1)=P(i+1)+epsilon;
        T(i,i+1)= -6.33E-03*Thalf(i,i+1,reservoir, fluid,numerical,P);           % Lower diagonal
        F2 =  T(i,i+1)*(P(i+1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        J(i,i+1) = (F2-F(i))/epsilon;
        P=Ptemp;
        
        J(i,i) = - (J(i,i+1)) + (B(i,i)/numerical.dt);   % Main diagonal
    elseif i==numerical.N
        T(i,i-1) = -6.33E-03*Thalf(i,i-1, reservoir, fluid,numerical,P);  
        F(i,1) = T(i,i-1)*(P(i-1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        
        Ptemp=P;
        P(i-1)=P(i-1)+epsilon;
        T(i,i-1)= -6.33E-03*Thalf(i,i-1,reservoir, fluid,numerical,P);            % Upper diagonal
        F2 = T(i,i-1)*(P(i-1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        J(i,i-1) = (F2-F(i))/epsilon;
        P=Ptemp;
        
        J(i,i) = - (J(i,i-1)) + (B(i,i)/numerical.dt);
    else        
        T(i,i-1)= -6.33E-03*Thalf(i,i-1,reservoir, fluid,numerical,P);            % Upper diagonal
        T(i,i+1)= -6.33E-03*Thalf(i,i+1,reservoir, fluid,numerical,P);           % Lower diagonal
        F(i,1) = T(i,i-1)*(P(i-1)-P(i))+ T(i,i+1)*(P(i+1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        
        Ptemp=P;
        P(i-1)=P(i-1)+epsilon;
        Ttemp= -6.33E-03*Thalf(i,i-1,reservoir, fluid,numerical,P);            % Upper diagonal
        F2 = Ttemp*(P(i-1)-P(i))+ T(i,i+1)*(P(i+1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        J(i,i-1) = (F2-F(i))/epsilon;
        P=Ptemp;
        
        Ptemp=P;
        P(i+1)=P(i+1)+epsilon;
        Ttemp= -6.33E-03*Thalf(i,i+1,reservoir, fluid,numerical,P);           % Lower diagonal
        F2 = T(i,i-1)*(P(i-1)-P(i))+ Ttemp*(P(i+1)-P(i)) + (B(i,i)/numerical.dt)*(P(i)-P_old(i)) - Q(i);
        J(i,i+1) = (F2-F(i))/epsilon;
        P=Ptemp;
        
        J(i,i) = - (J(i,i-1) + J(i,i+1)) + (B(i,i)/numerical.dt);   % Main diagonal
    end
   
end
%G= -0.433*6.33E-03*T*z;
%T=6.33E-03*T;