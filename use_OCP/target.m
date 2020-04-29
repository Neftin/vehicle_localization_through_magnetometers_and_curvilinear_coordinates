function [ value ] = target( z , auxdata )
% Evaluate the objective function

% Multiple shooting: evaluate function without integration

  N  = auxdata.N       ; % length of the simulation
  nx = auxdata.nx      ; % n-states
  np = auxdata.np      ; % n paramters (constant states)
  nu = auxdata.nu      ; % n-controls
  
  curvy_road = auxdata.curvy_road;
  
  magn = auxdata.magn; % magnetometers parameters
  
  NX = nx*N;
%  NP = np*N;
  NU = nu*(N-1);
      
  %% constants
  
  mi_0 = auxdata.mi_0;
  B_0  = auxdata.B_0;
  
  %% extract 
  
    x = reshape( z(1:NX) , [ nx , N ] ); % extract columns of states from z
    
    p_index = NX+1;             % first parameter
    
    p = z( p_index : p_index + np - 1 ); % parameters
    
    u_index = p_index + np; % first control
    
    u = reshape( z( u_index : u_index + NU - 1  ) , [ nu , N-1 ] ); % controls
    
    % convert in cartesian the curvilinear coordinates
    [ veh.x , veh.y , veh.psi ] = curvy_road.xypsi_by_snxi(  [ x(2,:) ;  x(3,:) ; x(4,:) ] );
      veh.m0                    = p(1:3);
      veh.D                     = p(4);
      veh.h                     = p(5);
  
  
  %% preallocate measures:

 
    for i = 1:magn.N

        meas(i).z      = zeros(length(veh.m0),N);        % measure preallocation

    end
    

  %% evaluate the measures

    for k=1:N

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure calculation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        rot         = makehgtform( 'zrotate', veh.psi(k) ); % True rotation of the vehicle in the ground
        m_k         = rot(1:3,1:3)*veh.m0;                  % Rotate the magnetic dipole (direct roation)
        D_k         = veh.D;                                % D fixed

        for j = 1:magn.N

            % k time index , j magnetometer index

            % construct r_k: relative position magnetometer -> dipole
            r_k            =  [ veh.x(k) - magn.x(j) , veh.y(k) - magn.y(j) , veh.h - magn.z(j) ]';

            rot_m          = makehgtform( 'zrotate', magn.th(j) ); % align magnetometer with the road

            % no noise measures
            meas(j).z(:,k) = MagnetoMeterSensor.MagnetometerSimulatedOutput( rot_m(1:3,1:3)'*r_k , rot_m(1:3,1:3)'*m_k , rot_m(1:3,1:3)'*B_0 , zeros(3)...
                , D_k )...
                ./mi_0;
            
            

        end


    end
    
  %% Error on the measures ( the target function )

    nm     = 3;
    
    errors = zeros( magn.N*nm , N );
    
    for i = 1:magn.N
        
        m_index = (i-1)*3 + 1;
        
        errors( m_index:(m_index + nm-1) , : ) = meas(i).z - auxdata.true_measures( m_index:(m_index + nm - 1) , : );

    end
    
    value = sum (errors(:).^2 ) ./ length( errors(:) );

end