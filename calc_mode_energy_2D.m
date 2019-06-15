function [intW, xInterp, yInterp] = calc_mode_energy_2D(w, x, y, r, minFeatureSize)
    % Finds mode energy of an antenna as a function of
    % of integration radius r. The antenna is assumed to have axial symmetry
    % thereby requiring only a 2D electric field profile.
    % w = 0.5 * Re( d[w*eps]/dw * |E|^2)
    % x, y are field coordinates (we are assuming Euclidean grid here)
    % r is assumed to be a vector with largest value occuring last
    % minFeatureSize is the min mesh discretization of the simulation.
    
    if r(end)>abs(x(end)) || r(end)>abs(y(end))
        error('Radius out of bounds')
    end

    intW = zeros(1,length(r));
    
    % Create interpolated points, to reduce error from defining circle,
    % and to reduce error from zero-ing energy near the antenna axis
    
    % Currently not automated, should be changed if necessary:
    interpSize1=1e-9; % feature size on interpolation mesh
    interpSize2=5e-11; % feature size in high density region (antenna)
    interpSize3=1e-11; % feature size right near axis
    
    xInterp=x(1);
    yInterp=y(1);
    
    for i=1:length(x)-1
        if abs(x(i)-x(i+1))>1.001*minFeatureSize
            xInterp=[xInterp;(x(i)+interpSize1:interpSize1:x(i+1)).']; 
        elseif x(i)==0 || x(i+1)==0
            xInterp=[xInterp;(x(i)+interpSize3:interpSize3:x(i+1)).'];
        else
            xInterp=[xInterp;(x(i)+interpSize2:interpSize2:x(i+1)).'];
        end
    end
    
    for j=1:length(y)-1
        if abs(y(j)-y(j+1))>1.001*minFeatureSize
            yInterp=[yInterp;(y(j)+interpSize1:interpSize1:y(j+1)).']; 
        elseif y(j)==0 || y(j+1)==0
            yInterp=[yInterp;(y(j)+interpSize3:interpSize3:y(j+1)).'];
        else
            yInterp=[yInterp;(y(j)+interpSize2:interpSize2:y(j+1)).'];
        end
    end
    
    [XInterp,YInterp] = meshgrid(xInterp,yInterp);
    
    [X,Y] = meshgrid(x,y);
    
    wInterp = interp2(X,Y,w.',XInterp,YInterp).';
    
    % integration goes as pi * int ( f(x,y) * x * dx * dy )
    % find values to multiply energy by
        
    scaleX = ones(length(xInterp),length(yInterp));

    for i=1:length(xInterp)
        scaleX(i,:)=abs(yInterp);
    end
    
    sphereFinder=(XInterp.^2+YInterp.^2).';
    
    for k=1:length(r)
        mask = sphereFinder <= r(k)^2;
        wScaled=mask.*scaleX.*wInterp;

        intW(1,k)=pi*trapz(xInterp,trapz(yInterp,wScaled,2));
        
    end
    
end
