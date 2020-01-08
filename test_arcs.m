
% test

ds = 1;
kv = linspace(0,1/3,1000);
theta_0 = 0;

a = curvilinear_arcs_trajectory( ds , theta_0 , kv );

s_rng = 0:ds/100:ds*length(kv);
[ x , y , th ] = a.xytheta_by_s( s_rng );

[ xv , yv , thv ] = a.xypsi_by_snxi( [ s_rng(1:10:end); 1*randn(1,length(s_rng(1:10:end))); zeros(1,length(s_rng(1:10:end))) ] );

a.fast_plot()
hold on;
plot(xv,yv,'redx');
hold off;
