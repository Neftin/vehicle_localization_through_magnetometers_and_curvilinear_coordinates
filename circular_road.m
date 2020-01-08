classdef circular_road
    %CIRCULAR_ROAD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R
    end
    
    methods
        function self = circular_road( radius )
            % CIRCULAR_ROAD Construct
            %   Detailed explanation goes here
            self.R = radius;
        end
        
        function theta = theta_from_s( self , s )
            % Return theta function of the curvilinear abscissa
            %   s has to be less than 2*pi*R and positive
            theta = (s ./ self.R + pi/2 );
            
        end
        
        function xy = xytheta_from_sn( s , n )
            
            x = (self.R-n)*cos( self.theta_from_s(s)-pi/2 );
            y = (self.R-n)*sin( self.theta_from_s(s)-pi/2 );

            xy = [ x y ]';
            
        end
        
    end
end

