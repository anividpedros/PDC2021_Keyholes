function C = PRV_2_DCM(t,e)

s = sin(t);
c = cos(t);
z = 1-c;

e1 = e(1);
e2 = e(2);
e3 = e(3);

C = [e1*e1*z+c e1*e2*z+e3*s e1*e3*z-e2*s;
    e2*e1*z-e3*s e2*e2*z+c e2*e3*z+e1*s;
    e3*e1*z + e2*s e3*e2*z-e1*s e3*e3*z+c];

end