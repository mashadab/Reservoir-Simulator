function [x,P] = Analytical_with_gravity (t,N,theta)
%input t:time (days),Scalar
%input N: number of grid blocks
%input theta: reservoir dip angle above the horizon in radians (positive..
%..for dipping up, negative for dipping down)
%output x: reservoir x-location (ft) 1*N vector 
%output P: pressure (psi) 1*N vector 
tol=eps;
rho_w=62.4;
L=4000; %Reservoir Length (ft)
za=2309; %Depth of left boundry of the reservoir (ft)
grid_vec=((1:N));
dx=L/N;
x=[(dx/2)*(2*grid_vec-1)];
phi = 0.25;
k = 15;
mu = 1;
cf = 10^-5;
dpg1=((za-L*sin(theta))*(rho_w/144))+14.6959502543;
P1 = 500-dpg1 ;
P0=0;
alpha = 6.33*10^-3*k/(phi*mu*cf);
n = 0;
SUM = 0;
for i=1:10000000
            new = (-1)^(n+1)*cos((n+1/2)*pi*x/L)*exp(-(2*n+1)^2/L^2*alpha*pi^2*t/4)/(2*n+1);
            SUM = 4/pi*(P1-P0).*new + SUM;
            error=max(abs(new));
            if error > tol
            n = n+1;
            else
                break
            end
end
dpg=((za-x*sin(theta))*(rho_w/144))+14.6959502543;
P= SUM +P1+dpg;
