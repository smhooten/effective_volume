function [VeffExt, Veff, VeffExtScaled, VeffScaled, PurcellF, wPeak] = calc_Veff(index, Ex, Ey, Ez, x, y, z, freq, r, minFeatureSize, peakEnergy, drudeModel, plots, Q)
    % Calculates effective mode volume given field data. Field data assumed
    % in the form Field(x,y,z,f)
    %
    % index is the complex refractive index, index(x,y,z,f)
    %
    % freq is the frequency (in Hz), can be a vector.
    %
    % r is a vector of radii to integrate over.
    %
    % Optional inputs:
    %
    % minFeatureSize is the minimum mesh size of the simulation.
    % This is only required if using a 2D field. (default=1)
    %
    % let w_e = 0.5 * Re( d(w*eps)/dw * |E|^2), energy density
    % peakEnergy is an integer flag that tells the program how to choose the location of
    % of peak energy density (i.e. Veff = total energy / peak energy density)
    % peakEnergy can be 0, 1, 2, or 3 (default=0)
    % 0: find w_e at center of simulation (most commonly used)
    % 1: find max(w_e) over entire simulation field data. This can be erroneous
    %    due to high field peaks near the surface of metal which are simulation artifacts
    % 2: find max(w_e) over entire simulation field data outside of metal regions. This can be 
    %    erroneous due to high field peaks near the surface of metal which are simulation artifacts
    % 3: find w_e averaged by 1 mesh cell at center of simulation. This probably should not
    %    not ever be used except if you suspect a source dipole was not properly subtracted.
    % 4: find max(w_e) over entire simulation, but the surface around metals is padded to
    %    avoid the problem mentioned above.
    % 5: User defined point to find peak w_e. Must edit the function.
    %
    % drudeModel is a boolean flag that tells the program whether to use a drude model to model the 
    % permittivity function of the metal. The drude model is typically too simplistic to get accurate results
    % (default=0)
    %
    % plots is a vector of integer flags that tells the program what to plot (default=[0])
    % 0: Plots Veff as a function of integration radius r, as well as y intercept extrapolation.
    %    All frequencies are plotted on the same curve
    % 1: Plots Veff as a function of frequency
    %
    % Q is the quality factor of the antenna for use with calculating the purcell factor given by
    % Fp = (3/(4*pi^2)) * Q / Veff * lambda^3
    % Note that if Q is not input, the Purcell factor calculation will be incorrect
    % (default=0)
    %
    % Outputs:
    %
    % VeffExt: the extracted mode volume in m^3 (vector with length of the input frequency vector)
    % Veff: the extracted mode volume in m^3 as a function of integration radius and input frequency
    % VeffExtScaled: the extracted mode volume normalized by lambda^3 (vector with length of the input frequency vector)
    % VeffScaled: the extracted mode volume normalized by lambda^3 as a function of integration radius and input frequency
    % PurcellF: Purcell factor (vector with length of the input frequency vector)
    % wPeak: extracted peak energy density (J/m^3), vector with length of the input frequency vector
    %
    % Note all units should be SI

    if ~exist('minFeatureSize','var')
        minFeatureSize=1;
    end
    
    if ~exist('peakEnergy','var')
        peakEnergy=0;
    end
    
    if ~exist('drudeModel','var')
        drudeModel=0;
    end

    if ~exist('plots','var')
        plots=[0];
    end
    
    if ~exist('Q','var')
        Q=0;
    end

    [~,eps0,~,c,~,~,~] = physical_constants();

    w=2*pi*freq;
    lambda=c./freq;

    % correct eps information for dispersive energy density

    [epsm, dwepsm] = permittivity_function(w, drudeModel);

    eps=index.^2;

    correctEps=abs((eps(:,:,:,1)-epsm(1))/epsm(1))<=0.1; % relative error less than 10%

    for i=1:length(freq)
        eps(:,:,:,i)=eps(:,:,:,i).*(~correctEps)+dwepsm(i)*correctEps; 

    end

    % find energy density

    E2=Ex.*conj(Ex)+Ey.*conj(Ey)+Ez.*conj(Ez);
    w_e=0.5*real(eps0*eps.*E2); 

    % find peak energy density

    if peakEnergy==0 % find max(w_e) at center of simulation
       wPeak=zeros(length(freq),1);
        
       xPeak=(length(x)+1)/2;
       yPeak=(length(y)+1)/2;
       zPeak=(length(z)+1)/2;
       yPeak = floor(yPeak)
       zPeak = floor(zPeak)
       wPeak(:,1)=w_e(xPeak,yPeak,zPeak,:);

    elseif peakEnergy==1 % find absolute max(w_e) over entire simulation
       wPeak=zeros(length(freq),1);

       for i=1:length(freq)
           wPeak(i,1)=max(max(max(w_e(:,:,:,i))));
       end

    elseif peakEnergy==2 % find absolute max(w_e) only in non-metallic regions, this should be used, but errors due to sharp features may be a problem
       w_e_Temp=zeros(size(w_e));
       wPeak=zeros(length(freq),1);

       for i=1:length(freq)
           w_e_Temp(:,:,:,i)=w_e(:,:,:,i).*(~correctEps);
           wPeak(i,1)=max(max(max(w_e_Temp(:,:,:,i))));
       end

    elseif peakEnergy==3 % find max(w_e) averaged near center
       xPeak=(length(x)+1)/2;
       yPeak=(length(y)+1)/2;
       zPeak=(length(z)+1)/2;

       wPeak=zeros(length(freq),1);
        
       if length(z)>1
           for i=1:length(freq)
               wPeak(i,1)=mean(mean(mean(w_e(xPeak-1:xPeak+1,yPeak-1:yPeak+1,zPeak-1:zPeak+1,i))));
           end
       else
           for i=1:length(freq)
               wPeak(i,1)=mean(mean(mean(w_e(xPeak-1:xPeak+1,yPeak-1:yPeak+1,1,i))));
           end
       end

    elseif peakEnergy==4 % find absolute max(w_e) only in non-metallic regions, where the boundary of the metal is fudged slightly to get rid of potentially untrue points
       w_e_Temp=zeros(size(w_e));
       wPeak=zeros(length(freq),1);
       correctEps2=correctEps;
       
       %how many points should be excluded?
       
       % two mesh cells
       if length(z)==1
           for i=3:length(x)-2
               for j=3:length(y)-2
                   if correctEps(i,j)
                       if any(any(~correctEps(i-2:i+2,j-2:j+2)))
                           correctEps2(i-2:i+2,j-2:j+2)=true;
                       end
                   end
               end
           end
       else
           for i=3:length(x)-2
               for j=3:length(y)-2
                   for k=3:length(z)-2
                       if correctEps(i,j,k)
                           if any(any(any(~correctEps(i-2:i+2,j-2:j+2,k-2:k+2))))
                               correctEps2(i-2:i+2,j-2:j+2,k-2:k+2)=true;
                           end
                       end
                   end
               end
           end
       end
       
       for i=1:length(freq)
           w_e_Temp(:,:,:,i)=w_e(:,:,:,i).*(~correctEps2);
           wPeak(i,1)=max(max(max(w_e_Temp(:,:,:,i))));
       end
       
    elseif peakEnergy==5
        %user specified peak
        for i=1:length(freq)
            wPeak(i,1)=w_e(379,301,301,i);
        end
        
    else
       error('not a valid peakEnergy flag')
    end


    % find total energy in simulation 
    
    if length(z)==1
        % 2D simulation (axial symmetry assumed)
        intW=zeros(length(freq),length(r));
        
        for i=1:length(freq)
            [intW(i,:),~,~]=calc_mode_energy_2D(w_e(:,:,:,i),x,y,r,minFeatureSize);
        end
        
    else
        % 3D
        intW=zeros(length(freq),length(r));
        
        for i=1:length(freq)
            intW(i,:)=calc_mode_energy_3D(w_e(:,:,:,i),x,y,z,r);
        end
    end
    
    % find effective volume
    
    Veff=zeros(length(freq),length(r));
    
    for i=1:length(r)
        Veff(:,i)=intW(:,i)./wPeak;
    end
    
    % extrapolate effective volume by subtracting y-intercept of linear
    % curve. NOTE this assumes that the integration radius was large enough
    % such that the antenna mode is completely accounted for. Note, there
    % also must be enough r data points such that the slope can be
    % calculated properly
    
    VeffExt=zeros(length(freq),1);
    
    for i=1:length(freq)
        m=(Veff(i,end)-Veff(i,end-1))/(r(end)-r(end-1));
        VeffExt(i)=Veff(i,end)-m*r(end);
        VeffExtScaled(i,:)=VeffExt(i)/lambda(i)^3;
    end
    
    % Calculate Purcell factor, given Q, assumes dipole resides on the axis
    % (x)
    
    PurcellF=zeros(length(freq),1);
    VeffScaled=zeros(length(freq),length(r));
    
    for i=1:length(freq)
        PurcellF(i)=(3/(4*pi^2))*Q/VeffExt(i)*lambda(i)^3;
        VeffScaled(i,:)=Veff(i,:)/lambda(i)^3;
    end

    % create plots
    
    for i=1:length(plots) 
        if plots(i)==0 % Veff as function of integration radius
            figure;
            
            for j=1:length(freq)
                plot(r*1e9,VeffScaled(j,:));
                hold on
            end
            title('Scaled Veff as Function of Integration Radius')
            xlabel('Integration Radius (nm)')
            ylabel('Veff / \lambda^3')
            hold off
            
            
            figure;

            plot(r*1e9,Veff/1e-18);
            title('Veff as Function of Integration Radius')
            xlabel('Integration Radius (nm)')
            ylabel('Veff (um^3)')

        elseif plots(i)==1 % Extrapolated Veff
            figure;
            
            plot(lambda*1e9,VeffExtScaled);
            title('Extrapolated Scaled Veff')
            xlabel('Wavelength (nm)')
            ylabel('Veff / \lambda^3')
            hold off
            
            figure;
            
            plot(lambda*1e9,VeffExt/1e-18)
            title('Extrapolated Veff')
            xlabel('Wavelength (nm)')
            ylabel('Veff (um^3)')
            
        end

    end

end
