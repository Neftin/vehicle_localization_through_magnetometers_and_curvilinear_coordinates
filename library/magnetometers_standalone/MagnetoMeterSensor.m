
% Class for magnetometers simulation

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

classdef MagnetoMeterSensor < handle
	% Model and notation are consistent with "extended target model" from the paper:
	% "Magnetometer Modeling  and validation for tracking metallic targets"

	% static version: no instantiation of the sensor needed

	methods( Static )
        
        % ~~~ Methods for generic sensors ~~~

		function [ yk ] = MagnetometerSimulatedOutput( r_k , m_k , B_0 , covM , D_k )
			% This method returns the output of the NOISY sensor (ONE READ ONLY)
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
			if any( size(r_k) ~= size(m_k) ) || (size(r_k,1) > 3) 
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

			% Noisy case
			if nargin >= 4 % if covM passed as input 
				if ~all(size(covM) == [3 3])
					error('Magnetometer:: wrong size of covariance matrix (must be 3x3)')
				else
					% generate noise and add (Gaussian noise)
					yk = yk + mvnrnd(zeros(1,3),covM)';
				end
            end
            
            % Noisy case
			if nargin >= 5 % if D passed as input
				if ( size(D_k,1) ~= 1 )
					error('Magnetometer:: wrong size of D must be row vector for n vehicles')
				else
					% add the rotation term
                    for i = length(D_k)
                        yk = yk + D_k(i)*B_0;
                    end
				end
			end
			

        end

% I don't remember the purpose of this function...

%         function [ yk ] = MicroMagnetometerSimulatedOutput(  r_k , m_k , B_0 , covM )
%             covM = zeros(3);
%             if nargin == 4 % if covM passed as input 
% 				if ~all(size(covM) == [3 3])
% 					error('Magnetometer:: wrong size of covariance matrix (must be 3x3)')
% 				else
% 					% generate noise and add (Gaussian noise)
% 					yk = yk + mvnrnd(zeros(1,3),covM)';
% 				end
% 			end
%             yk = MagnetoMeterSensor.MagnetometerSimulatedOutput( r_k , m_k , B_0 , covM )*10^6;
%         end
            

		% ~~~ Methods for specific sensors ~~~

		function [ yk ] = MMC5883MAoutput( r_k , m_k , B_0 , D )
			% Simulation of reading the sensor memsic MMC5883MA:
			% with data based on the datasheet

			% Note that the value provided is the RMS, which means
			% we have to extract the estimation of sigma^2: which is
			% rms_value * sqrt(n/n-1) = standard_deviation, but we use RMS since
			% we expect it was evaluated with n -> inf.

			covM = (12*10^-8)^2;
            
            if nargin == 4
                yk   = MagnetoMeterSensor.MagnetometerSimulatedOutput( r_k , m_k , B_0 , covM * eye(3) , D );
            else
                yk   = MagnetoMeterSensor.MagnetometerSimulatedOutput( r_k , m_k , B_0 , covM * eye(3) );
            end

			% quantization: 
			% range = 8 Gauss = 8x10^-4 Tesla
			% bit = 16, 2x8x10^-4 / 2^16 = 2.4414e-8

			quanto = 2.4414e-8;
			yk     = round(yk./quanto,0) .* quanto; % simulate quantization

        end
        
        function [ yk ] = MicroMMC5883MAoutput(  r_k , m_k , B_0 , D )
            if nargin == 4
                yk   = MagnetoMeterSensor.MMC5883MAoutput( r_k , m_k , B_0 , D );
            else
                yk   = MagnetoMeterSensor.MMC5883MAoutput( r_k , m_k , B_0  );
            end
        end
        



	end

end