function [ x , x_pred , P , P_pred , export ] = UKF_AN_step( fstate , x , P , hmeas , z , Q , R , params )
    % UKF   Unscented Kalman Filter for nonlinear dynamic systems with additive
    % noise

    % [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P 
    % for nonlinear dynamic system (for simplicity, noises are assumed as additive):

    %           x_k+1 = f(x_k) + w_k
    %           z_k   = h(x_k) + v_k
    % where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
    %       v ~ N(0,R) meaning v is gaussian noise with covariance R
    % Inputs:   f: function handle for f(x)
    %           x: "a priori" state estimate
    %           P: "a priori" estimated state covariance
    %           h: fanction handle for h(x)
    %           z: current measurement
    %           Q: process noise covariance 
    %           R: measurement noise covariance
    % Output:   x: "a posteriori" state estimate
    %           P: "a posteriori" state covariance
    %

    %% Computation of the basic parameters ( tunables and dimensions )
    % and weights

    L = numel(x);                               % numer of states

    m = numel(z);                               % numer of measurements

    alpha = 0.001;                              % default, tunable

    ki = 0;                                    % default, tunable

    beta = 2;                                   % default, tunable (2 optimum for Guassian)

    lambda = alpha^2*(L+ki)-L;                  % scaling factor

    c = L + lambda;                             % scaling factor

    %% computation of the weights ( depends on L and tunebales )
    
    W0    = lambda/c;
    W0c   = lambda/c + (1-alpha^2+beta);

    Wm    = [ W0  , ( 0.5/c )*ones(1,2*L) ];       % weights for means sum=1

    Wc    = [ W0c , ( 0.5/c )*ones(1,2*L) ];

    c = sqrt(c);

    %% Points

    X = sigmas(x,P,c);                            % sigma points around x
    
    export.X = X;                                 % only to plot

    %% unscented for f state

    [ x1 , X1 , P1 , X2 ] = ut( fstate , X , Wm , Wc , L , Q , params );       % unscented transformation of process

    % X1=sigmas(x1,P1,c);                         % sigma points around x1
    % X2=X1-x1(:,ones(1,size(X1,2)));             % deviation of X1
    
    x_pred = X1(:,1);
    P_pred = P1;

    [ z1 , Z1 , P2 , Z2 ] = ut( hmeas , X1 , Wm , Wc , m , R , params );       % unscented transformation of measurments

    P12 = X2 * diag(Wc) * Z2';                          % transformed cross-covariance

    K   = P12/P2;
    
    export.zRes  = (z-z1);
    export.zPred = z1; 
    export.K     = K;
    export.xRes  = K*(z-z1);
    
    x   = x1 + K*(z-z1);                                % state update

    P   = P1-K*P12';                                    % covariance update

function [ y , Y , P , Y1 ] = ut( f , X , Wm , Wc , n , R , params )

    %Unscented Transformation
    %Input:
    %        f: nonlinear map
    %        X: sigma points
    %       Wm: weights for mean
    %       Wc: weights for covraiance
    %        n: numer of outputs of f
    %        R: additive covariance
    %Output:
    %        y: transformed mean
    %        Y: transformed smapling points
    %        P: transformed covariance
    %       Y1: transformed deviations

    L = size(X,2);
    
    y = zeros(n,1);
    
    Y = zeros(n,L);
    
    % for each point apply the nonlinear map (model/measure) and computed
    % the mean
    for k = 1:L                   

        Y(:,k) = f( X(:,k) , params );  % In case of the model it is X+ the output, for measure is Y

        %y = y + Wm(k)*Y(:,k);           % transformed mean -> sum(Wm) = 1

    end
    
    y = Y(:,1); % !!! TOGLI 


    Y1 = Y - y(:)*ones(1,L);          % create vector L times the vector y

    P  = Y1 * diag(Wc) * Y1' + R;     % diag( wc(1)Y1(1)^2 + R(1) ... da 1..3 % Covariance of the prediction Pk-
    
function X = sigmas(x,P,c)
    % Sigma points around reference point
    % Inputs:
    %       x: reference point
    %       P: covariance
    %       c: coefficient sqrt(lambda+N)
    % Output:
    %       X: Sigma points
    
    X0 = x(:); % The mean

    DD_matrix = c * chol(P);
    
    XX_matrix = X0*ones(1,length(x)); 
    
    X = [ X0 , XX_matrix+DD_matrix , XX_matrix - DD_matrix ];







