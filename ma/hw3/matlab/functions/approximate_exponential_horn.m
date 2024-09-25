function [Z] = approximate_exponential_horn(l_exp_horn, a_0, m, delta, k, Z_l, rho, c)

S_tot = 1;

n_sections = l_exp_horn/delta;
for j = n_sections:-1:1
    r_1 = a_0*exp(m*delta*(j-1));
    r_2 = a_0*exp(m*delta*j);
    x_1 = delta*(r_1/(r_2-r_1));
    x_2 = x_1+delta;
    theta_1 = atan(k.*x_1);
    theta_2 = atan(k.*x_2);
    S_1 = r_1^2*pi;
    S_2 = r_2^2*pi;
    S_tot = S_1*S_tot;
    Z_n = 1i.*Z_l.*(sin(k.*delta-theta_2)./sin(theta_2))+(rho*c/S_2).*sin(k.*delta);
    Z_d = Z_l.*(sin(k.*delta+theta_1-theta_2)./(sin(theta_1).*sin(theta_2))) ...
        -(1i*rho*c/S_2).*(sin(k.*delta+theta_1)./sin(theta_1));
    Z = (rho*c/S_1).*(Z_n./Z_d);
    Z_l = Z;
end

end


