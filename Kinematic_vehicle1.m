function [ X ] = Kinematic_vehicle1( Xprec , params )

% State space discrete model for the Kinematic vehicle 
X  = zeros(9,1);

u  = Xprec(1);
s  = Xprec(2);
n  = Xprec(3);
xi = Xprec(4);

k  = params.k;    % instantaneous k
Ts = params.ts;   % sampling time


X(1)   = u;              % + noise here
X(2)   = s + Ts*u;
X(3)   = Ts*xi*u + n;
% various model yaw rate
X(4)   = xi;          % yaw rate = +ku -> il veicolo ha yaw rate esatto per rimanere allineato alla strada
%X(4)   = xi- Ts*k*u; % free flowing -> il veicolo va "per campi"

X(5:7) = Xprec(5:7);        % + noise
X(8)   = Xprec(8);          % + noise
X(9)   = Xprec(9);          % + noise

end

