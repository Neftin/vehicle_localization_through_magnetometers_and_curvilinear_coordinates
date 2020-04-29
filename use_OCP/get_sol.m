function [ sol_x , sol_p , sol_u ] = get_sol( z , auxdata )
% Parsing the z-vector to get the solution

    % usage:
    % [sol_x,sol_p,sol_u] = get_sol(sol_ocp,auxdata)
    % sol_x returns the state solution as column vectors
    % sol_p returns the parameters solution as column vectors
    % sol_u returns the control solution as column vectors

    N  = auxdata.N ; % length of the simulation
    nx = auxdata.nx; % n-states
    np = auxdata.np; % n paramters (constant states)
    nu = auxdata.nu; % n-controls

    sol_x         = reshape( z( 1 : nx*N ) , [ nx , N ] );
    z( 1 : nx*N ) = [];
    
    sol_p(:)      = z( 1 : np );
    sol_p         = sol_p(:);
    
    z( 1 : np )   = [];
    
    sol_u         = reshape( z( 1 : nu*( N - 1 ) ) , [ nu , ( N - 1 ) ] );


end