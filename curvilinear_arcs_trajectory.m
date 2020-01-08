classdef curvilinear_arcs_trajectory
    % curvilinear_arcs_trajectory
    % It store a trajectory with contants delta_S, and curvature constant for each step
    
    properties

        ds
        k
        theta 
        xy

    end
    
    methods
        function self = curvilinear_arcs_trajectory( ds , theta_zero , k_vector )

            if( any(size(ds) ~= [1,1]) )
                error('ds must be scalar')
            end

            self.ds       = ds;
            self.k        = k_vector;
            self.theta    = [ theta_zero , ( ds*cumsum( k_vector ) + theta_zero ) ];  % integral of piecewise constant curvature
            self.xy       = zeros(2,length(k_vector)); % cartesian coordinates preallocate

            for i = 2:length(k_vector)+1

                k = k_vector(i-1);

                self.xy(1,i) = self.xy(1,i-1) + (  - sin(self.theta(i-1)) + sin( k * ds + self.theta(i-1) ) ) / k ;
                if isnan(self.xy(1,i))
                    self.xy(1,i) = self.xy(1,i-1) + ds*cos(self.theta(i-1)); % if NaN means singularities -> linearize
                end

                self.xy(2,i) = self.xy(2,i-1) - (  - cos(self.theta(i-1)) + cos( k * ds + self.theta(i-1) ) ) / k ;
                if isnan(self.xy(2,i))
                    self.xy(2,i) = self.xy(2,i-1) + ds*sin(self.theta(i-1)); % if NaN means singularities -> linearize
                end

            end

        end

        function [ x,y,th ] = xytheta_by_s( self , s )

            s = s( s <= self.ds*length(self.k) ); % remove outbound values

            if any( mod(s,self.ds) ) % if any interstizial s (most of the cases)
                x  = zeros(1,length(s));
                y  = zeros(1,length(s));
                th = zeros(1,length(s));

                for i = 1:length(s)
                    idx = floor(s(i)./self.ds)+1; % get the right slice
                    dS  = mod(s(i),self.ds);      % get the extra portion of s
                    if dS > 0
                        kk    = self.k(idx);            % get curvature
                        x(i)  = self.xy(1,idx) + (  - sin(self.theta(idx)) + sin( kk * dS + self.theta(idx) ) ) / kk ; % integral on the arc
                        if isnan(x(i))
                            x(i) = self.xy(1,idx) + dS*cos(self.theta(idx)); % if NaN means singularities -> linearize
                        end
                        y(i)  = self.xy(2,idx) - (  - cos(self.theta(idx)) + cos( kk * dS + self.theta(idx) ) ) / kk ; % integral on the arc
                        if isnan(y(i))
                            y(i) = self.xy(2,idx) + dS*sin(self.theta(idx)); % if NaN means singularities -> linearize
                        end
                        th(i) = self.theta(idx) + dS*kk;

                    elseif dS == 0 % sharp on the junction among arcs
                        x(i)  = self.xy(1,idx);
                        y(i)  = self.xy(2,idx);
                        th(i) = self.theta(idx);

                    else
                        error('this simply cannot happen, write an email to giammarco(at)rocketmail.com')
                    end

                end

            else
                idx = (s./self.ds)+1
                x   = self.xy(1,idx);
                y   = self.xy(2,idx);
                th  = self.theta(idx);
            end

        end


        function [ x,y,psi ] = xypsi_by_snxi( self , curvicoord );
            % xi here is intended to be a generic angle with the tangent:
            % xi = psi-theta
            % s  = curvicoord(1,:);
            % n  = curvicoord(2,:);
            % xi = curvicoord(3,:);

            s  = curvicoord(1,:);
            n  = curvicoord(2,:);
            xi = curvicoord(3,:);

            s = s( s <= self.ds*length(self.k) ); % remove outbound values

            if any( mod(s,self.ds) ) % if any interstizial s (most of the cases)
                x   = zeros(1,length(s));
                y   = zeros(1,length(s));
                psi = zeros(1,length(s));
                th  = zeros(1,length(s));

                for i = 1:length(s)
                    idx = floor(s(i)./self.ds)+1; % get the right slice
                    dS  = mod(s(i),self.ds);      % get the extra portion of s
                    if dS > 0
                        kk     = self.k(idx);            % get curvature

                        th(i)  = self.theta(idx) + dS*kk;

                        psi(i) = th(i) + xi(i); % get Psi

                        x(i)  = self.xy(1,idx) + (  - sin(self.theta(idx)) + sin( kk * dS + self.theta(idx) ) ) / kk ; % integral on the arc
                        if isnan(x(i))
                            x(i) = self.xy(1,idx) + dS*cos(self.theta(idx)); % if NaN means singularities -> linearize
                        end
                        x(i) = x(i) - n(i)*sin( th(i) ); % add n


                        y(i)  = self.xy(2,idx) - (  - cos(self.theta(idx)) + cos( kk * dS + self.theta(idx) ) ) / kk ; % integral on the arc
                        if isnan(y(i))
                            y(i) = self.xy(2,idx) + dS*sin(self.theta(idx)); % if NaN means singularities -> linearize
                        end
                        y(i) = y(i) + n(i)*cos( th(i) ); % add n

                    elseif dS == 0 % sharp on the junction among arcs
                        th(i)  = self.theta(idx);
                        psi(i) = th(i) + xi(i);
                        x(i)   = self.xy(1,idx) - n(i).*sin( th(i) );
                        y(i)   = self.xy(2,idx) + n(i).*cos( th(i) );
                    else
                        error('this simply cannot happen, write an email to giammarco(at)rocketmail.com')
                    end

                end

            else
                idx = (s./self.ds)+1;
                th  = self.theta(idx);
                psi = th + xi;
                x   = self.xy(1,idx) - n.*sin(th);
                y   = self.xy(2,idx) + n.*cos(th);
                
            end

        end
        
        function [ k_out  ] = get_k_from_s( self , s )
            % return the curvature at s
            idx    = floor(s./self.ds)+1;
            k_out  = self.k( idx );
        end

        function fast_plot( self )
            s_rng = 0:self.ds/100:self.ds*length( self.k );
            [ x , y , th ] = self.xytheta_by_s( s_rng );
            plot(x,y,'blacko');
            xlabel('x (m)');
            ylabel('y (m)');
            axis equal;
        end
        
    end
    
end