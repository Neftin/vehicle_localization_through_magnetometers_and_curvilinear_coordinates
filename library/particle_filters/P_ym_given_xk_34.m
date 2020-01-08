function P_ym = P_ym_given_xk_34( yk , x_propagated , sigma_m , muM , Km )
    % function to generate the particles from the sensor reading.
    % Based on article: 
    % Vehicle Tracking Based on Fusion of Magnetometer and Accelerometer 
    %         Sensor Measurements With Particle Filtering
    
    % offset of the earth magnetic field has to be removed
    
    if ( any(size(yk) ~= [ 3 1 ]) ) || ( any(size( x_propagated ) ~= [3 1]) ) 
        error('x and yk in "P_ym_given_xk_34()" must be 3x1 vectors')
    end

    r       = x_propagated;      % distance

    Hm      = (10^-7)*( 3 * (r * r.') - norm(r)^2 * eye(3) ) ./ norm(r)^5; 
    % Sensor measure function (nonlinear), according to eq(34) 

    Km      = Km;                                   % COVARIANCE of noise: obtainable by Monte Carlo
    mu      = Hm * muM;                             % PRIOR of the measurement given the prior of the  Dipole moment
    covM    = sigma_m.^2 .* ( eye(3) + Hm*Km*Hm.'); % Resulting gaussian co-variance
    P_ym    = mvnpdf( yk , mu , covM );             % it's x to span, not yk: drawn from P(Yk|Xk)
   
    
end
