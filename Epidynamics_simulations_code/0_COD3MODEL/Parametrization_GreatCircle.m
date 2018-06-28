function [mu2, mu1, nu]=Parametrization_GreatCircle(E,F,R,z)

mu2=R*(E(1)*cos(z)+F(1)*sin(z));
mu1=-R*(E(2)*cos(z)+F(2)*sin(z));
nu=R*(E(3)*cos(z)+F(3)*sin(z));

end