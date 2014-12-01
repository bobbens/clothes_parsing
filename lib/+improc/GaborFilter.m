classdef GaborFilter
    %GABORFILTER
    % Naive implementation of Gabor filter
    % Code taken from http://en.wikipedia.org/wiki/Gabor_filter
    % Recursive filter should be implemented to save computational cost
    %
    % Usage:
    %   [ g ] = GaborFilter(...)
    %   [ Y ] = g.apply( X )
    %
    % Kota Yamaguchi
    
    properties
        sigma  = 1.41        % Standard deviation
        theta  = (0:3)/4*pi  % Orientation
        lambda = 16          % Wavelength
        psi    = 0.5*pi      % Phase
        gamma  = 1           % Aspect ratio
    end
    
    methods
        function [ this ] = GaborFilter( varargin )
            %GABORFILTER
            props = properties(this);
            for i = 1:2:length(varargin)
                ind = strcmp(varargin{i},props);
                if any(ind), this.(props{ind}) = varargin{i+1}; end
            end
        end
        
        function [ Y ] = apply( this, X )
            %APPLY
            if size(X,3)>1, X = im2double(rgb2gray(X)); end
            Y = zeros(size(X,1),size(X,2),this.size);
            F = this.filters();
            for i = 1:this.size
                Y(:,:,i) = filter2(F{i},X,'same');
            end
        end
        
        function [ Y ] = filters( this )
            %FILTERS return cell array of filter banks
            Y = cell(1,this.size);
            i = 1;
            for s = this.sigma
                for t = this.theta
                    for l = this.lambda
                        for p = this.psi
                            for g = this.gamma
                                Y{i} = this.create_filter(s,t,l,p,g);
                                i = i+1;
                            end
                        end
                    end
                end
            end
        end
        
        function [ Y ] = size( this )
            %SIZE
            Y = length(this.sigma) * length(this.theta) *...
                length(this.lambda) * length(this.psi) *...
                length(this.gamma);
        end
    end
    
    methods(Static)
        function gb = create_filter(sigma,theta,lambda,psi,gamma)
            %GABOR_FUNC
            sigma_x = sigma;
            sigma_y = sigma/gamma;

            % Bounding box
            nstds = 3;
            xmax = max(abs(nstds*sigma_x*cos(theta)),...
                       abs(nstds*sigma_y*sin(theta)));
            xmax = ceil(max(1,xmax));
            ymax = max(abs(nstds*sigma_x*sin(theta)),...
                       abs(nstds*sigma_y*cos(theta)));
            ymax = ceil(max(1,ymax));
            xmin = -xmax; ymin = -ymax;
            [x,y] = meshgrid(xmin:xmax,ymin:ymax);
            % Rotation 
            x_theta= x*cos(theta)+y*sin(theta);
            y_theta=-x*sin(theta)+y*cos(theta);

            gb= 1/(2*pi*sigma_x *sigma_y) *...
                exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)) .*...
                cos(2*pi/lambda*x_theta+psi);
        end
    end
    
end

