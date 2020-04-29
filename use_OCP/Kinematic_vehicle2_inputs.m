function [ X ] = Kinematic_vehicle2_inputs( Xprec , Uprec , auxdata )

% State space discrete model for the Kinematic vehicle 
X  = zeros(auxdata.nx,1);

u       = Xprec(1);
s       = Xprec(2);
n       = Xprec(3);
xi      = Xprec(4);
yawrate = Xprec(5); % yawraterate input

Ts = auxdata.Ts;   % sampling time

X(1)   = u + Uprec(1)*Ts; % acc input
X(2)   = s + Ts*u;
X(3)   = Ts*xi*u + n;
% various model yaw rate:
X(4)   = Ts*yawrate + xi;   % yaw rate = +ku -> il veicolo ha yaw rate esatto per rimanere allineato alla strada (+ rumore)
% X(4)   = xi - Ts*k*u; % free flowing -> il veicolo va "per campi"
X(5)   = yawrate + Uprec(2)*Ts; % yaw rate + yawraterate inpus

end