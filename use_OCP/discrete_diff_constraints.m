function h_vector = discrete_diff_constraints( z , auxdata , f )
% return the constraints as a function of the variables:
% discrete systems case x+ = f(x,u)

% Multiple shooting: evaluate function without integration

  N  = auxdata.N ; % length of the simulation
  nx = auxdata.nx; % n-states
  np = auxdata.np; % n paramters (constant states)
  nu = auxdata.nu; % n-controls
  ts = auxdata.Ts;       % sampling time

  NX = nx*N;
%  NP = np*N;
  NU = nu*(N-1);
  
  h_vector = zeros(NX-nx,1); % constraint vector ( h(x) )
  % h(x) = 0 means h(x) =  x^+ - f(x)
  
  for i = 2:N
      
      
      rng_x      = (i-2)*nx+1 : (i-1)*nx; % x
      rng_x_plus = (i-1)*nx+1 : (i)*nx  ; % x+ 
                    
      index_u = nx*N+np;  
      rng_u   = index_u + ( i-2 )*nu + (1:nu);
      
      x      = z( rng_x      );
      x_plus = z( rng_x_plus );
      u      = z( rng_u );
      
      h_vector( rng_x ) = x_plus - f( x , u , auxdata );   % function of the system x+ = f(x,u)
      
  end
  
%         fprintf('Norm of constraint evaluated: %f\n', norm(h_vector(1)))
  
end
