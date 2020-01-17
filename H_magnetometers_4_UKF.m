function z_meas  = H_magnetometers_4_UKF( X , params )
    % model of the measurement function: multiple sensor model
    % it works for only one set of states: x column vector
    
    N      = length(params.s0); % how many sensors?
    z_meas = zeros(3*N,1);      % vector of measurements
    
    mi_0 = 4*pi*1e-07; 
    
    % Unpack parameters (vectors: multiple sensors)
    s0  = params.s0;
    n0  = params.n0;
    z0  = params.z0;
    B0  = params.B0./mi_0;  % B0 has to be IN THE GROUND REFERENCE FRAME -> and normaòlized on mi_0
    th0 = params.th0;
    
   
    % Unpack the state
    if length(X) == 10 % model 1
        s     = X(2); % second state: curvilinear abscissa
        n     = X(3); % third state:  curvilinear ordinata
        xi    = X(4); % fourth state: curvilinear heading

        z     = X(10); % !!! check the number of the state it can be a source of errors
        m_veh = X(6:8);
        D     = X(9);  % !!! check the number of the state it can be a source of errors
    elseif length(X) == 9
        s     = X(2); % second state: curvilinear abscissa
        n     = X(3); % third state:  curvilinear ordinata
        xi    = X(4); % fourth state: curvilinear heading

        z     = X(9); % !!! check the number of the state it can be a source of errors
        m_veh = X(5:7);
        D     = X(8); % !!! check the number of the state it can be a source of errors
    else
        error('unmanaged dimension of the state X in the UFK')
    end
    
    % Curvilinear rotation
    rot   = makehgtform('zrotate',xi);
    rot   = rot(1:3,1:3);
   
    % Iterate over sensors
    for i = 1:N
        
        r_k = [ s - s0(i) , n - n0(i) , z - z0(i) ]';
        m_k = rot*m_veh(:); % rotate vehicle only  ( curvilinear coordinates already contains theta )
        D_k = D;
        
        rotG = makehgtform('zrotate',th0(i)); % !!! check the rotation direction
        rotG = rotG(1:3,1:3);
        
        
    
        B_0  = rotG'*B0(:); % inverse rotation for the earth magnetic field
        
        idx    = ( (i-1)*3 + 1 ) : ( i*3 ); % 1..3 4..6 and so on
        z_meas(idx) = H_magnetometer( r_k , m_k , B_0 , D_k );
        
    end

end


function yk      = H_magnetometer( r_k , m_k , B_0 , D_k )
    % This method returns the output of the ideal sensor
    % given the set of the dipole moments around it and the earth magnetic field 
    % note: in order to do multiple reading a loop in needed for now (it iterates only
    % through dipoles of targets)
    %
    % Usage:
    % r_k  = [ r_x(:) ; r_y(:) ; r_z(:) ] % column vectors (concatenate for multiple dipoles)
    % m_k  = [ m_x(:) ; m_y(:) ; m_z(:) ] % column vectors (concatenate for multiple dipoles)
    % B_0  = [ B_0_x  ; B_0_y  ; B_0_z  ] % column vector
    % D_k  = [ D ]                        % scalar for fixed field
    % (D*B0)
    % covM = []                           % covariance matrix (3x3)

    % Put check on arguments
    if any( size(r_k) ~= size(m_k) ) || (size(r_k,1) > 3 ) 
        error('MagnetoMeter:: wrong size of some argument')
    end
    if size(r_k) == [ 0 0 ]
        r_k = zeros(3,1);
        m_k = zeros(3,1);
    end


    % Number of dipoles and space dimension
    d = size(r_k,2);
    n = size(r_k,1);

    % compute the matrices JM of each dipole

    mi_0 = 4*pi*1e-07; 
    
    mi_0 = 1; % NORMALIZE !!!

    for i = 1:d
    % computation for each 
        rk = r_k(:,i);
        mk = m_k(:,i);
        I  = eye(n);

        % equation 2c
        Jm(:,:,i) = ( mi_0 / ( (4*pi)*( norm(rk)^5 ) ) )  *  ( 3*( rk * rk' ) - ( rk' * rk ) * I );

        % compute fields
        Jmmk(:,i) = Jm(:,:,i)*mk;

    end

    yk = B_0 + sum(Jmmk,2); % sum over the rows to get the total moment 

    if nargin >= 4 % if D passed as input
        if ( size(D_k,1) ~= 1 )
            error('Magnetometer:: wrong size of D. It must be row vector for n vehicles')
        else
            % add the rotation term
            for i = length(D_k)
                yk = yk + D_k(i)*B_0;
            end
        end
    end


end