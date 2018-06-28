%   Model in Saggio et al. 2017 JMN
%   [x1dot;x2dot;zdot]=[x,y,z];
%   b = 1       here (for different values of b cfr Dumortier et al. 91)
%   c -->       inverse of timescale separation constant or velocity of slow variable
%   dstar -->   determines when z has to change direction
%   R -->       radius in unfolding. Use R=0.4 by default
%   E, F -->    vectors describing the arc of great circle (the path on
%   which the system moves)
%   N -->       N=1 means the upper branch of the equilibrium manifold is
%   the resting state. This holds for most of the classes. N=2 for the
%   lower branch


function xdot = Cod3_Model(~,x,~,~,c,R,dstar,E,F,N)

    % for a given value of z (x(3)), gives position along the path
    [mu2, mu1, nu]=Parametrization_GreatCircle(E,F,R,x(3));

    % for a given position on the map, gives values of resting (or silent) state 
    x_rs=real(eval_resting_state_cartesian(mu2,mu1,nu,N));

    % System
    x1dot = - x(2);
    x2dot = -( -x(1)^3 +mu2*x(1) +mu1 + x(2)*( nu + x(1) + x(1)^2));
    zdot =  -c*(sqrt((x(1)-x_rs)^2+x(2)^2)-dstar);
    
xdot = [x1dot;x2dot;zdot];
