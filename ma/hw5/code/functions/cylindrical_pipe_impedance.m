function [Z] = cylindrical_pipe_impedance(l, k, r, rho, c)

S = pi.*r.^2;
Z_0 = rho*c./S;
Z = 1i*Z_0*tan(k.*l);

end