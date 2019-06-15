function [hbar, eps0, q, c, mu0, Z0, euler] = physical_constants()
    % Useful physical constants

    hbar = 6.63e-34 / (2*pi);
    eps0 = 8.854e-12;
    q = 1.602e-19;
    c = 299792458;
    mu0 = 4*pi*1e-7;
    Z0 = sqrt(mu0/eps0);
    euler = 0.577215664901532;

end
