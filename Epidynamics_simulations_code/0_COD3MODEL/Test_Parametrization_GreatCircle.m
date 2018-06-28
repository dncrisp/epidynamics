% By setting A and B allows to visualize the complete great circle in the
% unfolding

% A and B determine the path
A=[0.2731,-0.05494,0.287]; 
B=[0.3331,0.074,0.2087];

% general settings
R=0.4; % radius of the unfolding
z=0:0.01:(2*pi); % full great circle (2 pi)

% compute path's vectors
E=A/R;
F=cross(cross(A,B),A);
F=F/norm(F);

% compute path in cartesian coordinates
mu2=R*(E(1)*cos(z)+F(1)*sin(z));
mu1=-R*(E(2)*cos(z)+F(2)*sin(z));
nu=R*(E(3)*cos(z)+F(3)*sin(z));

% plot curve followed in parameter space
openfig('FOCUS_r0p4/Unfolding_r0p4.fig')
hold on
plot3(mu2,-mu1,nu,'r')
hold off
