
%% execute after curvytrack.m in the main folder ...

N = length(rng_flt);
nx = 5;
nu = 2;
np = 5;

auxdata.N        = N;      % number of cells or intervals
auxdata.nx       = nx;     % number of states
auxdata.nu       = nu;     % number of controls
%auxdata.nb       = nb;     % number of boundary conditions
auxdata.np       = np;
auxdata.magn     = magn;
auxdata.B_0      = B_0;
auxdata.mi_0     = mi_0;
auxdata.Ts       = ts;

  NX = nx*N;
%  NP = np*N;
  NU = nu*(N-1);

%% put the measures and the road in auxdata

auxdata.true_measures = zV;
auxdata.curvy_road = curvy_road;


%% Initial guess

x0_ocp = zeros(nx,N); 
u0_ocp = zeros(nu,N-1); 
p0_ocp = zeros(np,1);

extender = ones(1,N);

x0_ocp(1,:) = x_init( 1 , extender );
x0_ocp(2,:) = linspace( magn.s(1)-2,magn.s(end)+2,N );
x0_ocp(3,:) = 0;
x0_ocp(4,:) = 0;
x0_ocp(5,:) = 0;

p0_ocp(1)  = 80;
p0_ocp(2)  = 80;
p0_ocp(3)  = 80;
p0_ocp(4)  = 1;
p0_ocp(5) = 0.3;

u0_ocp(1,:) = 0;
u0_ocp(2,:) = 0;

z0_ocp = [ reshape( x0_ocp , [1 , NX] ) , p0_ocp(:)' , reshape( u0_ocp , [1 , NU] ) ]';

%% try target with the true states (should be zero)

z_example = xpu2z( xTrue(1:5,rng_flt)+[ 0 50 10 0 0 ]' , xTrue(6:10,1) , u0_ocp , auxdata );
error_example      = target( z_example , auxdata );
[ ~ , constraint_example ] = constraints( z_example , auxdata , @Kinematic_vehicle2_inputs );

%% !!! SE BLOCCA
% call NLP to solve discretized ocp  
options = optimoptions('fmincon','Display','iter-detailed','Algorithm','sqp','StepTolerance',1e-5,'MaxIterations',1000,'FiniteDifferenceType','central','DiffMinChange',1e-3,'DiffMaxChange',1,'MaxFunctionEvaluations',10^6);
     
%% Solve OCP

sol = fmincon( @(z)target(z,auxdata) , ...
        z0_ocp, ...
        [],[],[],[],[],[], ...
        @(z)constraints( z , auxdata , @Kinematic_vehicle2_inputs ) , ...
        options);

%% get solutions        
        
[sol_x,sol_p,sol_u] = get_sol( sol , auxdata );
 

%% plot solution

% convert in cartesian the curvilinear coordinates
[ from_ocp.x , from_ocp.y , from_ocp.psi ] = curvy_road.xypsi_by_snxi(  [ sol_x(2,:) ;  sol_x(3,:) ; sol_x(4,:) ] );
  from_ocp.m0                    = sol_p(1:3);
  from_ocp.D                     = sol_p(4);
  from_ocp.h                     = sol_p(5);
  
figure(3)
dr.sdraw = 0:0.001:road.len;
[ dr.xb0 , dr.yb0 , ~ ] = curvy_road.xytheta_by_s( dr.sdraw );
[ dr.xbl , dr.ybl , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ;  road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
[ dr.xbr , dr.ybr , ~ ] = curvy_road.xypsi_by_snxi( [ dr.sdraw ; -road.w*ones(1,length(dr.sdraw)) ; zeros(1,length(dr.sdraw)) ]);
plot(dr.xb0,dr.yb0,'black--',dr.xbl,dr.ybl,'blue-',dr.xbr,dr.ybr,'blue-');
axis equal;
grid on;
hold on;

figure(3)
plot(from_ocp.x,from_ocp.y,'red','LineWidth',2)
plot(veh.x,veh.y,'blue:','LineWidth',2);
plot(magn.x,magn.y,'blackx');
hold off;
      
      
%%

function z = xpu2z( x,p,u,auxdata )

  N  = auxdata.N       ; % length of the simulation
  nx = auxdata.nx      ; % n-states
  np = auxdata.np      ; % n paramters (constant states)
  nu = auxdata.nu      ; % n-controls
  NX = nx*N            ;
  NU = nu*(N-1)        ;
  
    z = [ reshape( x , [1 , NX] ) , p(:)' , reshape( u , [1 , NU] ) ]';
end

          