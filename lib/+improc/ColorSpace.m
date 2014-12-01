classdef ColorSpace
    %COLORSPACE color conversion functions
    % Code taken from http://www.easyrgb.com/index.php
    % Kota Yamaguchi 2011
    
    methods(Static)
        function [ Y ] = xyz2rgb( X )
            %XYZ2RGB
            X = X / 100;
            Y = cat(3,...
                3.2406*X(:,:,1) - 1.5372*X(:,:,2) - 0.4986*X(:,:,3),...
               -0.9689*X(:,:,1) + 1.8758*X(:,:,2) + 0.0415*X(:,:,3),...
                0.0557*X(:,:,1) - 0.2040*X(:,:,2) + 1.0570*X(:,:,3));
            for i = 1:3
                y = Y(:,:,i);
                ind = y > 0.0031308;
                y(ind) = 1.055 * (y(ind).^(1/2.4)) - 0.055;
                y(~ind) = 12.92 * y(~ind);
                Y(:,:,i) = y;
            end
            Y = im2uint8(Y);
        end
        
        function [ Y ] = rgb2xyz( X )
            %RGB2XYZ
            X = im2double(X);
            for i = 1:3
                x = X(:,:,i);
                ind = x > 0.04045;
                x(ind) = ((x(ind) + 0.055)/1.055).^2.4;
                x(~ind) = x(~ind) / 12.92;
                X(:,:,i) = x * 100;
            end
            Y = cat(3,...
                0.4124*X(:,:,1) + 0.3576*X(:,:,2) + 0.1805*X(:,:,3),...
                0.2126*X(:,:,1) + 0.7152*X(:,:,2) + 0.0722*X(:,:,3),...
                0.0193*X(:,:,1) + 0.1192*X(:,:,2) + 0.9505*X(:,:,3));
        end
        
        function [ Y ] = xyz2lab( X, ref )
            %XYZ2LAB from XYZ to CIE-L*ab
            if nargin < 2, ref = [95.047,100.000,108.883]; end
            for i = 1:3
                x = X(:,:,i) ./ ref(i);
                ind = x > 0.008856;
                x(ind) = x(ind).^(1/3);
                x(~ind) = 7.787*x(~ind) + (16/116);
                X(:,:,i) = x;
            end
            Y = cat(3,...
                116*X(:,:,2) - 16,...
                500*(X(:,:,1)-X(:,:,2)),...
                200*(X(:,:,2)-X(:,:,3)));
        end
        
        function [ Y ] = lab2xyz( X, ref )
            %LAB2XYZ from CIE-L*ab to XYZ
            if nargin < 2, ref = [95.047,100.000,108.883]; end
            y = (X(:,:,1)+16) / 116;
            Y = cat(3, X(:,:,2)/500 + y, y, y - X(:,:,3)/200);
            for i = 1:3
                y = Y(:,:,i);
                y_3 = y.^3;
                ind = y_3 > 0.008856;
                y(ind) = y_3(ind);
                y(~ind) = (y(~ind)-16/116) / 7.787;
                Y(:,:,i) = y * ref(i);
            end
        end
        
        function [ Y ] = rgb2lab( X, ref )
            %RGB2LAB from RGB to CIE-L*ab
            if nargin < 2, ref = [95.047,100.000,108.883]; end
            Y = improc.ColorSpace.rgb2xyz(X);
            Y = improc.ColorSpace.xyz2lab(Y, ref);
        end
        
        function [ Y ] = lab2rgb( X, ref )
            %LAB2RGB from CIE-L*ab to RGB
            if nargin < 2, ref = [95.047,100.000,108.883]; end
            Y = improc.ColorSpace.lab2xyz(X, ref);
            Y = improc.ColorSpace.xyz2rgb(Y);
        end
        
        function [ Y ] = rgb2hsv( X )
            %RGB2HSV
            X = im2double(X);
            [r,g,b] = deal(X(:,:,1),X(:,:,2),X(:,:,3));
            [d_r,d_g,d_b] = deal(zeros(size(X,1),size(X,2)));
            [H,S] = deal(zeros(size(X,1),size(X,2)));
            v_min = min(X,[],3);
            v_max = max(X,[],3);
            d_max = v_max - v_min;
            V = v_max;
            ind = d_max ~= 0;
            d_max = d_max(ind);
            v_max = v_max(ind);
            S(ind) = d_max ./ v_max;
            d_r(ind) = (((v_max-r(ind))./6)+(d_max/2))./d_max;
            d_g(ind) = (((v_max-g(ind))./6)+(d_max/2))./d_max;
            d_b(ind) = (((v_max-b(ind))./6)+(d_max/2))./d_max;
            I = r==V&ind; H(I) =       d_b(I)-d_g(I);
            I = g==V&ind; H(I) = (1/3)+d_r(I)-d_b(I);
            I = b==V&ind; H(I) = (2/3)+d_g(I)-d_r(I);
            H(H<0) = H(H<0)+1;
            H(H>1) = H(H>1)-1;
            Y = cat(3,H,S,V);
        end
        
        function [ Y ] = hsv2rgb( X )
            %HSV2RGB
            [H,S,V] = deal(X(:,:,1),X(:,:,2),X(:,:,3));
            [R,G,B] = deal(zeros(size(X,1),size(X,2)));
            ind = S~=0;
            R(~ind) = V(~ind); G(~ind) = V(~ind); B(~ind) = V(~ind);
            h = H * 6; h(h==6) = 0; i = floor(h);
            v1 = V.*(1-S);
            v2 = V.*(1-S.*(h-i));
            v3 = V.*(1-S.*(1-(h-i)));
            I = i==0&ind; R(I) =  V(I); G(I) = v3(I); B(I) = v1(I);
            I = i==1&ind; R(I) = v2(I); G(I) =  V(I); B(I) = v1(I);
            I = i==2&ind; R(I) = v1(I); G(I) =  V(I); B(I) = v3(I);
            I = i==3&ind; R(I) = v1(I); G(I) = v2(I); B(I) =  V(I);
            I = i==4&ind; R(I) = v3(I); G(I) = v1(I); B(I) =  V(I);
            I = i==5&ind; R(I) =  V(I); G(I) = v1(I); B(I) = v2(I);
            Y = im2uint8(cat(3,R,G,B));
        end
    end
end