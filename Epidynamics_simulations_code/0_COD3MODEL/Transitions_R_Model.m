

function xdot = Transitions_R_Model(~,x,~,~,c,R,E,F,dstar,c2,N)

Ru=R+x(4);

mu2=Ru*(E(1)*cos(x(3))+F(1)*sin(x(3)));
mu1=-Ru*(E(2)*cos(x(3))+F(2)*sin(x(3)));
nu=Ru*(E(3)*cos(x(3))+F(3)*sin(x(3)));

    x_rs=real(eval_resting_state_cartesian(mu2,mu1,nu,N));
   
    x1dot = - x(2);
    x2dot = -( -x(1)^3 +mu2*x(1) +mu1 + x(2)*( nu + x(1) + x(1)^2));
    
    zdot =  -c*(sqrt((x(1)-x_rs)^2+x(2)^2)-dstar);

    udot= c2;
   
    xdot = [x1dot;x2dot;zdot;udot];
end
