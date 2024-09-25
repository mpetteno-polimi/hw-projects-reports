function [Z] = conical_horn_impedance(l, k, r_t, r_m, Z_l, rho, c)

x_1 = l.*(r_t./(r_m-r_t));
x_2 = x_1+l;
theta_1 = atan(k.*x_1);
theta_2 = atan(k.*x_2);
S_1 = r_t.^2.*pi;
S_2 = r_m.^2.*pi;
Z_n = 1i.*Z_l.*(sin(k.*l-theta_2)./sin(theta_2))+(rho*c./S_2).*sin(k.*l);
Z_d = Z_l.*(sin(k.*l+theta_1-theta_2)./(sin(theta_1).*sin(theta_2))) ...
    -(1i*rho*c./S_2).*(sin(k.*l+theta_1)./sin(theta_1));
Z = (rho*c./S_1).*(Z_n./Z_d);

end