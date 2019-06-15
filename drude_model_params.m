function [wp, tau, epsInf] = drudeModelParams(material)
    % Provides wp, tau, and epsInf for use in calculating
    % the permittivity function of a metal.

    % Parameters here are based on fitted lumerical model
    
    [hbar,~,q,~,~,~,~] = physicalConstants;

    if strcmp(material, 'Silver') || strcmp(material, 'Ag')
        wp = 8.9*q/hbar;
        tau = 15e-15;
        epsInf = 37 + 2*(1i);
    
    else
        disp('material not supported');

    end
