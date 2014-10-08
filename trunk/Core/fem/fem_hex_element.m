classdef fem_hex_element < fem_element
    %FEM_HEX_ELEMENT Linear hexahedral element
    %
    % Node indices are assumed to be ordered as follows:
    %   1:4, nodes of one face ordered clockwise w.r.t. outward normal
    %   5:8, corresponding nodes on opposite face (will be CCW)
    %
    %   Copyright 2013 C. Antonio Sanchez [antonios@ece.ubc.ca]
    
     %% Constants
    properties (Constant, Access=private)
                
        INTEGRATION_COORDS_GAUSS_8 = ...
           [[ -1, -1, -1;
               1, -1, -1;
               1,  1, -1;
              -1,  1, -1;
              -1, -1,  1;
               1, -1,  1;
               1,  1,  1;
              -1,  1,  1 ]/sqrt(3) ones(8,1) ];
        
        
        % cubature point locations
        ipntLength = 8;     % number of cubature points
        % cubature weights
        ipntW = fem_hex_element.INTEGRATION_COORDS_GAUSS_8(:,4)';  
        % cubature point locations
        ipntLoc = fem_hex_element.INTEGRATION_COORDS_GAUSS_8(:,1:3);    
        
        % number of nodes in hex elements
        hexNumNodes = 8;
    end
    
    %% Private properties
    properties (Access = private)
        ipntRules = {};
    end
    
    %% Static methods
    methods (Static)
        function [pnts, w] = getCubatureRule( order )
            % Computes a Gaussian cubature rule over the hex element
            % 
            % [pnts, w] = getCubatureRule(order)
            %   Computes and returns the Gaussian cubature points and
            %   weights to accurately integrate a function of the given
            %   order
            
            [p, a] = fem_element.getGaussQuadrature(order);
            n = length(p);
            N = n^3;
            pnts = zeros(N,3);
            w = zeros(1, N);
            
            idx = 1;
            for i=1:n
                for j=1:n
                    for k=1:n
                        pnts(idx,:) = [locs(i), locs(j), locs(k)];
                        w(idx) = a(i)*a(j)*a(k);
                        idx = idx+1;
                    end
                end
            end
            
        end
        
        function N = getShapeFunction(xi, eta, mu)
            % Static shape function routine
            %
            % N = getShapeFunction(xi, eta, mu)
            %   Evaluates the element's shape function at natural
            %   coordintes [xi, eta, mu]
            
            xm = 1-xi;
            xp = 1+xi;
            em = 1-eta;
            ep = 1+eta;
            mm = 1-mu;
            mp = 1+mu;
            
            N = [xm.*em.*mm, xp.*em.*mm, xp.*ep.*mm, xm.*ep.*mm, ...
                 xm.*em.*mp, xp.*em.*mp, xp.*ep.*mp, xm.*ep.*mp]/8;
        end
    end
    
    %% Implemented methods
    methods
        function n = getNumNodes(~)
            n = fem_hex_element.hexNumNodes;
        end
        
        function N = getN(~,xi,eta,mu)
            xm = 1-xi;
            xp = 1+xi;
            em = 1-eta;
            ep = 1+eta;
            mm = 1-mu;
            mp = 1+mu;
            
            N = [xm.*em.*mm, xp.*em.*mm, xp.*ep.*mm, xm.*ep.*mm, ...
                 xm.*em.*mp, xp.*em.*mp, xp.*ep.*mp, xm.*ep.*mp]/8;
        end
        
        function dN = getdN(~,xi,eta,mu)
            xm = 1-xi;
            xp = 1+xi;
            em = 1-eta;
            ep = 1+eta;
            mm = 1-mu;
            mp = 1+mu;
            
            dN = [-em*mm, em*mm, ep*mm,-ep*mm,-em*mp, em*mp, ep*mp,-ep*mp;
                  -xm*mm,-xp*mm, xp*mm, xm*mm,-xm*mp,-xp*mp, xp*mp, xm*mp;
                  -xm*em,-xp*em,-xp*ep,-xm*ep, xm*em, xp*em, xp*ep, xm*ep];
            dN = dN/8;
        end
        
        function [in] = isNaturalInside(~, xem)
            in = (xem(1)>=-1 && xem(1)<=1 && xem(2)>=-1 && xem(2)<=1 ...
                    && xem(3)>=-1 && xem(3)<=1);
        end
        
        function [ipnts, w] = getIPnts(~)
            ipnts = fem_hex_element.ipntLoc;
            w = fem_hex_element.ipntW;
        end
        
    end
    
    %% Over-ridden methods
    methods
         function B = getB(~, dNdx)
            B = [dNdx(1,1) 0 0  dNdx(1,2) 0 0  dNdx(1,3) 0 0 ...
                    dNdx(1,4) 0 0  dNdx(1,5) 0 0  dNdx(1,6) 0 0 ...
                    dNdx(1,7) 0 0  dNdx(1,8) 0 0;
                 0 dNdx(2,1) 0  0 dNdx(2,2) 0  0 dNdx(2,3) 0 ...
                    0 dNdx(2,4) 0  0 dNdx(2,5) 0  0 dNdx(2,6) 0 ...
                    0 dNdx(2,7) 0  0 dNdx(2,8) 0;
                 0 0 dNdx(3,1)  0 0 dNdx(3,2)  0 0 dNdx(3,3) ...
                    0 0 dNdx(3,4)  0 0 dNdx(3,5)  0 0 dNdx(3,6) ...
                    0 0 dNdx(3,7)  0 0 dNdx(3,8);
                 dNdx(2,1) dNdx(1,1) 0  dNdx(2,2) dNdx(1,2) 0 ...
                    dNdx(2,3) dNdx(1,3) 0  dNdx(2,4) dNdx(1,4) 0 ...
                    dNdx(2,5) dNdx(1,5) 0  dNdx(2,6) dNdx(1,6) 0 ...
                    dNdx(2,7) dNdx(1,7) 0  dNdx(2,8) dNdx(1,8) 0;
                 dNdx(3,1) 0 dNdx(1,1)  dNdx(3,2) 0 dNdx(1,2) ...
                    dNdx(3,3) 0 dNdx(1,3)  dNdx(3,4) 0 dNdx(1,4) ...
                    dNdx(3,5) 0 dNdx(1,5)  dNdx(3,6) 0 dNdx(1,6) ...
                    dNdx(3,7) 0 dNdx(1,7)  dNdx(3,8) 0 dNdx(1,8);
                 0 dNdx(3,1) dNdx(2,1)  0 dNdx(3,2) dNdx(2,2) ...
                    0 dNdx(3,3) dNdx(2,3)  0 dNdx(3,4) dNdx(2,4) ...
                    0 dNdx(3,5) dNdx(2,5)  0 dNdx(3,6) dNdx(2,6) ...
                    0 dNdx(3,7) dNdx(2,7)  0 dNdx(3,8) dNdx(2,8)];
        end
    end
    
    %% Hex-specific methods
    methods
        function hex = fem_hex_element(E)
        % Constructor, creates a single or array of elements
        %
        % hex = fem_hex_element(E)
        %   Creates a set of M hexahedral elements given an Mx8 matrix of
        %   node indices, arranged CW w.r.t. the outward normal
        
            if (nargin == 0) 
                hex.nodeIdxs = [];
            else
                if (size(E,1)<2)
                    hex.nodeIdxs = E;
                else
                    % pre-allocate array
                    hex(size(E,1)) = fem_hex_element([]);
                    for i=1:size(E,1)
                        hex(i).nodeIdxs = E(i,:);
                    end
                end
            end
        end
    end
    
end

