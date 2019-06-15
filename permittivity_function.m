function [epsm, dwepsm] = permittivity_function(w, drudeModel, material)
    % Provides the permittivity, epsm, and the relevant derivative
    % for use in finding the energy density in a metal at frequency w
    % dwepsm = d[w*epsm]/dw

    if drudeModel % Don't use, inaccurate
        [wp, tau, epsInf] = drude_model_params(material);
        
        epsm = epsInf-wp^2./(w.^2+(1i)*w/tau);
        dwepsm = epsInf+wp^2*(w.^2-1/tau^2)./((w.^2+1/tau^2).^2);

	% Alternate number:
        % dwepsm = real(epsm)+2*w.*imag(epsm)*tau; % Koenderick/Ruppin paper
        
    else % Palik data for silver at 1550nm
        epsm = -87+(1i)*8.74;
        dwepsm = 94-(1i)*11.5;

    end

end
