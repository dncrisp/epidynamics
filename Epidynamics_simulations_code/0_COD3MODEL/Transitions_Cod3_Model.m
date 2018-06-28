% Marisa Saggio 2016

% model for transitions between bursting classes (Saggio et al. 2017, JMN)
%   [x1dot;x2dot;zdot]=Cod3_Model(~,x,~,~,c,R,dstar,E,F,N)
%   b = 1       here (for different values of b cfr Dumortier et al. 91)
%   c -->       inverse of timescale separation constant or velocity of slow variable
%   dstar -->   determines when z has to change direction
%   R -->       radius in unfolding. Use R=0.4 by default
%   G,H,L,M -->    vectors describing the path on
%   which the system moves)
%   cA, cB -->  velocities at which Ain moves towards Aend and Bin towards
%   Bend
%   N -->       N=1 means the upper branch of the equilibrium manifold is
%   the resting state. This holds for most of the classes. N=2 for the
%   lower branch

function xdot = Transitions_Cod3_Model(~,x,~,~,c,R,G,H,L,M,dstar,cA,cB,N)

    % A1 goes to B1 because of z (following great circle)
    % A1 goes to A2 and B1 goes to B2 because of uA and uB (following great circle)

    Au=R*(G*cos(x(4))+H*sin(x(4))); % A(u) and B(u)
    Bu=R*(L*cos(x(5))+M*sin(x(5)));
    
    Eu=Au/R;                        % vectors for path from A(u) to B(u)
    Fu=cross(cross(Au,Bu),Au);
    Fu=Fu/norm(Fu);
   
    mu2=R*(Eu(1)*cos(x(3))+Fu(1)*sin(x(3))); 
    mu1=-R*(Eu(2)*cos(x(3))+Fu(2)*sin(x(3)));
    nu=R*(Eu(3)*cos(x(3))+Fu(3)*sin(x(3)));

    x_rs=real(eval_resting_state_cartesian(mu2,mu1,nu,N));

    % System
   
    x1dot = - x(2);
    x2dot = -( -x(1)^3 +mu2*x(1) +mu1 + x(2)*( nu + x(1) + x(1)^2));
    zdot =  -c*(sqrt((x(1)-x_rs)^2+x(2)^2)-dstar);
    uAdot= cA;
    uBdot= cB;
   
    xdot = [x1dot;x2dot;zdot;uAdot;uBdot];
end
