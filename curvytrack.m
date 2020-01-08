%% Curvylinear Kalman for magnetometers
% function definitions

%  #######################################################
%  #  _______   __  _____ ______ _______ _ ____   _____  #
%  #         \ |  |   ___|   ___|__   __| |    \ |       #
%  #          \|  |  __| |  __|    | |  | |     \|       #
%  #       |\     | |____| |       | |  | |   |\         #
%  #  _____| \____| _____|_|       |_|  |_|___| \______  #
%  #                                                     #
%  #######################################################

% Coded by Giammarco Valenti
% /// KEEP IT SIMPLE! \\\

% reference papers:
% [1] Single target tracking using vector magnetometers
% [2] Magnetometer Modeling and Validation for Tracking Metallic Targets
    % for magnetometer modelling
% [3] Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
% Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
% pp.401-422, 2004. 


close all;
clear all; 
addpath(genpath('../library'));
rng('default')
% states:
% [ u , s , n , xi , ms , mn , mz , D , z ]
% u  = longitudinal velocity no slip
% s  = curvilinear abscissa
% n  = normal so s
% xi = angle respect to tanget of the road
% ms,n,z = dipole moment vector (one for now)
% D constant of the vehicle [3] for rotaion magnetic
% z height

%% Road geometry definition

road.len          = 100; % m
road.w            = 2;   % m (half width)
road.ds           = 1;
road.min_R        = 50;
road.k            = [ linspace(0,1./road.min_R,road.len/(2*road.ds)) linspace(1./road.min_R,-1./road.min_R,road.len/(2*road.ds)) ];
theta_0           = 0;

curvy_road = curvilinear_arcs_trajectory( road.ds , theta_0 , road.k );

%% plot road
% 
figure(1)
dr.sdraw = 0:0.001:road.len;
[ dr.xb0 , dr.yb0 , ~ ] = curvy_road.xytheta_by_s( dr.sdraw );
[ dr.xbl , dr.ybl , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ;  road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
[ dr.xbr , dr.ybr , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ; -road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
plot(dr.xb0,dr.yb0,'black--',dr.xbl,dr.ybl,'blue-',dr.xbr,dr.ybr,'blue-');
axis equal;
grid on;
hold on;

clear dr

%% Create vehicle data

ts       = 0.005; % 200 Hz
veh.Tsin = 7;
veh.v    = 10;
veh.s    = 0;
veh.xi   = 0;
veh.n    = 0;
t        = 0;

while( veh.s(end) < road.len )
    veh.v(end+1)   = veh.v(1) + 2*sin( 1*pi*(t(end)/veh.Tsin) );
    veh.s(end+1)   = veh.v(end)*ts + veh.s(end);
    veh.n(end+1)   = 1 + 0.5*sin(2*pi*(t(end)/veh.Tsin));
    veh.xi(end+1)  = atan2( veh.n(end) - veh.n(end-1) , veh.s(end) - veh.s(end-1) ); % !!! VERIFICA
    t(end+1)       = t(end)+ts;
end

% remove outbounds (time stops one value before)
veh.v(end)   = [];
veh.s(end)   = [];
veh.n(end)   = [];
veh.xi(end)  = [];
t(end)       = [];



% figure(1)
% convert in cartesian
[ veh.x , veh.y , veh.psi ] = curvy_road.xypsi_by_snxi(  [ veh.s ;  veh.n ; veh.xi ] );
% plot(veh.x,veh.y,'red');

veh.z        = 0.4*ones( size(veh.x) );

%% Magnetometer positions

magn.s      = 30:4:46; % it determines also the number of couples of magnetometers 
magn.line_n = [2 -2];
magn.n      = magn.line_n'*ones(1,length(magn.s)); % use it as magn(1,:) etc
magn.xi     = zeros(1,length(magn.s));
magn.N      = length(magn.s)*length(magn.line_n); % total magnetometers



for i = 1:length(magn.line_n)
    % fill line by line from 1 to N
    
    idx = i:magn.line_n:magn.N;
    
    [ magn.x( idx ) , magn.y( idx ) , magn.th( idx ) ] = curvy_road.xypsi_by_snxi( [magn.s ; magn.n(i,:) ; magn.xi ] );
    
end

magn.z      = zeros(size(magn.x));

figure(1);

for i = 1:magn.N
    
    plot(magn.x(i),magn.y(i),'blackx');
    hold on;
    text(magn.x(i)+1,magn.y(i),num2str(i)); % number of the magnetometer
    
end

%% Magnetic quantities

B_0 = [ 17 17 2 ]'*1e-06;   % Tesla

veh.D  = 1.01;               % reasonable value
veh.m0 = [ -200 100 -200 ]'; % for now one dipole only

warehouse.m_ground = zeros(length(veh.m0),length(t));                % for now only one dipole !!! add dimension
warehouse.r_k      = zeros(length(magn.N),length(veh.m0),length(t)); % for now only one dipole !!! add dimension  

% preallocate measures
for i = 1:magn.N
    
    meas(i).z      = zeros(length(veh.m0),length(t));        % measure preallocation
    
end

    meas_test.z    = zeros(length(veh.m0)*magn.N,length(t)); % To test the model (and the handle of the function)

%% simulate measures
% compute measures -> simulation
for k=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % measure calculation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rot         = makehgtform( 'zrotate', veh.psi(k) ); % True rotation of the vehicle in the ground
    m_k         = rot(1:3,1:3)*veh.m0;                  % Rotate the magnetic dipole (direct roation)
    D_k         = veh.D;                                % D fixed
    
    warehouse.m_ground(:,k) = m_k;

    for j = 1:magn.N
        
        % k time index , j magnetometer index
        
        % construct r_k: relative position magnetometer -> dipole
        r_k            =  [ veh.x(k) - magn.x(j) , veh.y(k) - magn.y(j) , veh.z(i) - magn.z(j) ]';
        
        warehouse.r_k(j,:,k) = r_k; % relative positions
        
        rot_m          = makehgtform( 'zrotate', magn.th(j) ); % align magnetometer with the road
        
        meas(j).z(:,k) = MagnetoMeterSensor.MMC5883MAoutput( rot_m(1:3,1:3)'*r_k , rot_m(1:3,1:3)'*m_k , rot_m(1:3,1:3)'*B_0 , D_k );  
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    % TEST function of the multiple sensor
    %%%%%%%%%%%%%%%%%%%%%%
    
    X_true    = [ veh.v(k) , veh.s(k) , veh.n(k) , veh.xi(k) , veh.m0' , veh.D , veh.z(k) ]';
    
    % this step is important has to be formalized !!!
    tmp        = [ magn.s ; magn.s ];
    params.s0  = tmp(:)';            % two by two the same s
    params.n0  = magn.n(:)';         % two by two DIFFERENT n
    params.z0  = zeros(1,magn.N);
    params.th0 = magn.th(:)'; % this come from the knoeledge of the curvilinear coordinates !!! (it can be estimated from B_0)
    params.B0  = B_0;
    
    meas_test.z(:,k) = H_magnetometers_4_UKF( X_true , params ); % test the function used in UKF

end

%% plot measures/test

figure(2)

for i=1:magn.N
    
    plot( veh.s , meas(i).z );
    hold on;
    idx_test = ((i-1)*3+1) : (i*3);
    plot( veh.s , meas_test.z( idx_test , : ),'Color',[.3 .3 .3 .5],'linewidth',3);
    
end

%% filtering parameters

Qinit = diag( [ 5     , 5   , 2   , 0.01    , 50 , 50 , 50 , 0.1 , 0.1  ].^2 ); % Error on initial conditions
Q     = diag( [ 1*ts , 0.1 , 0.1 ,   0.5*ts  , 50 , 50 , 50 , 0.1 , 0.01 ].^2 ); % Q additive
R     = ( 12*10^-3 ).^2*eye(magn.N*3);

%% Filtering

% filtering only around the magnetometer region
rng_flt = find( veh.s > ( min(magn.s) -2 ) & veh.s < ( max(magn.s) +5 ) );


% DEBUG QUANTITIES !!!
sample_of_a_sampling    = zeros(2,9,length( rng_flt )); % contains diagonl elements of the samples
the_values_of_a_warrior = zeros(9,length( rng_flt ));   % eigenvalues of prediction covariance.

% initial states:
x_init = [ 10 , ( min(magn.s) -2 + 0.5*randn(1) ) , 1 , 0 ,  -150 , 80 , -150 , 1 , 0.3 ];

n_x = length(x_init);
n_z = magn.N*3;

x      = x_init; %initial state
P      = Qinit;
xV     = zeros(n_x,length(rng_flt));
sV     = zeros(n_x,length(rng_flt));
zV     = zeros(n_z,length(rng_flt));
zRes   = zeros(n_z,length(rng_flt));
xRes   = zeros(n_x,length(rng_flt));
zPred  = zeros(n_z,length(rng_flt));

disp('filtering')

for k=1:length(rng_flt)
   
    true_index = rng_flt(k); % index of the real time 

    % collect measurements from all the sensors
    z = zeros(3,magn.N); 
    for i=1:magn.N
      z(:,i) = meas(i).z( : , true_index );   
    end
    
    z = z(:);

    zV(:,k) = z;                             % save measurment
   
    % FIXED PARAMETERS
    tmp        = [ magn.s ; magn.s ];
    params.s0  = tmp(:)';
    params.n0  = magn.n(:)';
    params.z0  = zeros(1,magn.N);
    params.th0 = magn.th(:)'; % this come from the knowledge of the curvilinear coordinates !!! (in the reality it can be estimated from B_0)
    params.B0  = B_0;
    % CHANGING PARAMETERS
    params.k   = curvy_road.get_k_from_s( x(2) );
    params.ts  = ts;
    

    % UKF %%%%%%
    %%%%%%%%%%%%
    [x, xpred , P , Ppred , export ]  = UKF_AN_step( @Kinematic_vehicle1 , x , P , @H_magnetometers_4_UKF , z , Q , R , params );
    %%%%%%%%%%%%

    % DEBUG prediction only
    %     
    %     x = xpred;
    %     P = Ppred;
    %     
    %     sample_of_a_sampling(1,:,k)  = diag(export.X(1:9,2:10));
    %     sample_of_a_sampling(2,:,k)  = diag(export.X(1:9,11:19));
    %     the_values_of_a_warrior(:,k) = eig(Ppred);

    xV(:,k)     = x;                            % save estimate
    zRes(:,k)   = export.zRes;                  % save residuals
    kV{k}       = export.K;                     % save kalman gains
    zPred(:,k)  = export.zPred;                 % save zPred
    xRes(:,k)   = export.xRes;
  
end

disp('end filtering')

%% plot

figure(1)
dr.sdraw = 0:0.001:road.len;
[ dr.xb0 , dr.yb0 , ~ ] = curvy_road.xytheta_by_s( dr.sdraw );
[ dr.xbl , dr.ybl , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ;  road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
[ dr.xbr , dr.ybr , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ; -road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
plot(dr.xb0,dr.yb0,'black--',dr.xbl,dr.ybl,'blue-',dr.xbr,dr.ybr,'blue-');
axis equal;
grid on;
hold on;

figure(1)
decimation = 50;
[ final_traj.x , final_traj.y , final_traj.th ] = curvy_road.xypsi_by_snxi( xV(2:4,1:decimation:end) );
quiver(final_traj.x,final_traj.y,1*cos(final_traj.th),1*sin(final_traj.th),'red');
plot(veh.x(1:decimation:end),veh.y(1:decimation:end),'blue');
plot(magn.x,magn.y,'blackx');
hold off;

%axis([ (min(magn.x)-5) (max(magn.x)+5) (min(magn.y)-5) (max(magn.y)+5)  ]); 
grid minor;

%%
figure(3)
subplot(1,2,1)
plot( [ veh.xi(rng_flt) ; xV( 4 , : ) ]' );
legend('true xi','estimation xi')
subplot(1,2,2)
plot( [ veh.n(rng_flt) ; xV( 3 , : ) ]' );
legend('true n','estimation n')

grid minor;

%% Other elaboraions

% how the Kalman gain ( line of the angle ) evolves...

% columns the lines of xi

for i=1:length(rng_flt)
    columns_of_kV_of_xi(:,i) = kV{i}(:,4);
end

figure(4)
yyaxis right
plot( vecnorm(columns_of_kV_of_xi) );
set(gca, 'YScale', 'log')
yyaxis left
plot( [ veh.xi(rng_flt) ; xV( 4 , : ) ]' );
legend('true','estimate')

title('Norm of column vector of the xi-column of Kalman gain')
xlabel('samples')
hold off;




for i=1:length(rng_flt)
    columns_of_kV_of_xi(:,i) = kV{i}(:,3);
end

figure(5)
yyaxis right
plot( vecnorm(columns_of_kV_of_xi) );
set(gca, 'YScale', 'log')
yyaxis left
plot( [ veh.n(rng_flt) ; xV( 3 , : ) ]' );
legend('true','estimate')

title('Norm of column vector of the n-column of Kalman gain')
xlabel('samples')
hold off;










