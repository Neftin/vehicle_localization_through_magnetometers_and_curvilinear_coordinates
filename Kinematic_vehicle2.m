function [ X ] = Kinematic_vehicle2( Xprec , params )

% State space discrete model for the Kinematic vehicle 
X  = zeros(9,1);

u       = Xprec(1);
s       = Xprec(2);
n       = Xprec(3);
xi      = Xprec(4);
yawrate = Xprec(5);

k  = params.k;    % instantaneous k
Ts = params.ts;   % sampling time


X(1)   = u;              % + noise here
X(2)   = s + Ts*u;
X(3)   = Ts*xi*u + n;
% various model yaw rate:
X(4)   = Ts*yawrate + xi;   % yaw rate = +ku -> il veicolo ha yaw rate esatto per rimanere allineato alla strada (+ rumore)
% X(4)   = xi - Ts*k*u; % free flowing -> il veicolo va "per campi"
X(5)   = yawrate; % yaw rate + noise

X(6:8) = Xprec(6:8);        % + noise
X(9)   = Xprec(9);          % + noise
X(10)  = Xprec(10);         % + noise

end