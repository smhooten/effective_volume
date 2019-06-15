function intW = calc_mode_energy_2D(w, x, y, z, r)
    % Finds total energy of a mode with 3D energy density, w,
    % as a function of integration radius, r.
    % intW is approximately int(4*pi*r^2 w dr) over limits 0 to r
    % w = 0.5 * Re( d[w*eps]/dw * eps0 * |E|^2) is assumed.
    % x, y, z are vectors of field coordinates
    % r is assumed to be a vector with largest value occuring last
    
    if r(end)>abs(x(end)) || r(end)>abs(y(end)) || r(end)>abs(z(end))
        error('Radius out of bounds')
    end
    
    intW = zeros(1,length(r));

    [X,Y,Z] = meshgrid(x,y,z);
    
    sphereFinder = X.^2+Y.^2+Z.^2;
    sphereFinder = permute(sphereFinder,[2 1 3]);
    
    for k=1:length(r)
        mask = sphereFinder <= r(k)^2;
        wMask=mask.*w;

        intW(1,k)=trapz(x,trapz(y,trapz(z,wMask,3),2),1);
    end
    
end
