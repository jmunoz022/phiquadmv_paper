function [eps,beta,u0,x1,x2,y1,y2,T] = data(P)
%
% DATA - Function that returns the data of the problem
%
% OUTPUT: 
%   eps           - diffusion coefficient 
%   beta          - advection vector
%   u0            - initial condition
%   [x1 x2 y1 y2] - limits of the rectangular domain in space
%   T             - final time
%   P             - problem number:
%                    - P=1 for the Heat equation
%                    - P=2 for the Eriksson-Johnson probem
%                    - P=3 for the Hochburck-Ostermann equation
%                    - P=4 for the Allen-Cahn equation
%
if P == 1 %Heat equation
    eps = 1;
    beta = [0 0];
    u0 = @(x,y) sin(pi*x).*sin(pi*y);
    x1 = 0;
    x2 = 1;
    y1 = 0;
    y2 = 1;
    T = 1;
elseif P==2 %Eriksson-Johnson probem
    eps = 10^-2;
    beta = [1 0];
    C = 10;
    r1 = (1-sqrt(1+4*(pi*eps)^2))/(2*eps);
    r2 = (1+sqrt(1+4*(pi*eps)^2))/(2*eps);
    u0 = @(x,y) C*x.*(y.^2-0.25)+cos(pi*y).*(exp(r1*x)-exp(r2*x))./(exp(-r1)-exp(-r2));
    x1 = -1;
    x2 = 0;
    y1 = -0.5;
    y2 = 0.5;
    T = 1;
elseif P == 3 %Hochburck-Ostermann equation
    eps = 1;
    beta = [0 0];
    u0 = @(x,y) x.*(1-x).*y.*(1-y);
    x1 = 0;
    x2 = 1;
    y1 = 0;
    y2 = 1;
    T = 1;
elseif P == 4 %Allen-Cahn equation
    eps = 1;
    beta = [0 0];
    N = 7;
    alpha = 0.75;
    theta = @(x,y) atan((y-0.5)./(x-0.5)).*(x>0.5)+(pi+atan((y-0.5)./(x-0.5))).*(x<=0.5);
    u0 = @(x,y) tanh((0.25+0.1*cos(N*theta(x,y))-sqrt((x-0.5).^2+(y-0.5).^2))/(sqrt(2)*alpha));
    x1 = 0;
    x2 = 1;
    y1 = 0;
    y2 = 1;
    T = 2e-2;
end

end