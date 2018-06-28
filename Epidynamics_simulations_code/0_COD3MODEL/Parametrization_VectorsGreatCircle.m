function [E,F]=Parametrization_VectorsGreatCircle(p1,p2,R)

E=p1/R;
F=cross(cross(p1,p2),p1);
F=F/norm(F);

end