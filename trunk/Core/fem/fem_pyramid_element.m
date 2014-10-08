classdef fem_pyramid_element < fem_element
    %FEM_PYRAMID_ELEMENT Linear pyramid element
    %
    % Node indices are assumed to be ordered as follows:
    %   1:4, nodes of the square face
    %   5, the final node
    %
    %   Copyright 2013 C. Antonio Sanchez [antonios@ece.ubc.ca]
    
    %% Constants
    properties (Constant, Access=private)
        
        INTEGRATION_COORDS_GAUSS_5 = ...
           [[-1,-1;
              1,-1;
              1, 1;
             -1, 1] * (8/5)*sqrt(2/15), -ones(4,1)*2/3, ones(4,1)*81/100;
             0, 0, 2/5, 125/27];
        
        % cubature point locations
        ipntLength = 5;     % number of cubature points
        % cubature weights
        ipntW = fem_pyramid_element.INTEGRATION_COORDS_GAUSS_5(:,4)';  
        % cubature point locations
        ipntLoc = fem_pyramid_element.INTEGRATION_COORDS_GAUSS_5(:,1:3);    
        
        % number of nodes in pyramid elements
        pyramidNumNodes = 5;
    end
    
    %% Implemented methods
    methods
        function n = getNumNodes(~)
            n = fem_pyramid_element.pyramidNumNodes;
        end
        
        function N = getN(~,xi,eta,mu)
            xm = (1-xi);
            xp = (1+xi);
            em = (1-eta);
            ep = (1+eta);
            mm = (1-mu);
            mp = (1+mu);
            
            N = [xm.*em.*mm, xp.*em.*mm, xp.*ep.*mm, xm.*ep.*mm, 4*mp]/8;
        end
        
        function dN = getdN(~,xi,eta,mu)
            xm = (1-xi);
            xp = (1+xi);
            em = (1-eta);
            ep = (1+eta);
            mm = (1-mu);
            
            dN = [-em*mm, em*mm, ep*mm,-ep*mm, 0;
                  -xm*mm,-xp*mm, xp*mm, xm*mm, 0;
                  -xm*em,-xp*em,-xp*ep,-xm*ep, 4]/8;
        end
        
        function [in] = isNaturalInside(~, xem)
            in = (xem(1)>=-1 && xem(1)<=1 && xem(2)>=-1 && xem(2)<=1 ...
                    && xem(3)>=-1 && xem(3)<=1);
        end
        
        function [ipnts, w] = getIPnts(~)
            ipnts = fem_pyramid_element.ipntLoc;
            w = fem_pyramid_element.ipntW;
        end
        
    end
    
    %% Over-ridden methods
    methods
         function B = getB(~, dNdx)
            B = [dNdx(1,1) 0 0  dNdx(1,2) 0 0  dNdx(1,3) 0 0 ...
                    dNdx(1,4) 0 0  dNdx(1,5) 0 0;
                 0 dNdx(2,1) 0  0 dNdx(2,2) 0  0 dNdx(2,3) 0 ...
                    0 dNdx(2,4) 0  0 dNdx(2,5) 0;
                 0 0 dNdx(3,1)  0 0 dNdx(3,2)  0 0 dNdx(3,3) ...
                    0 0 dNdx(3,4)  0 0 dNdx(3,5);
                 dNdx(2,1) dNdx(1,1) 0  dNdx(2,2) dNdx(1,2) 0 ...
                    dNdx(2,3) dNdx(1,3) 0  dNdx(2,4) dNdx(1,4) 0 ...
                    dNdx(2,5) dNdx(1,5) 0;
                 dNdx(3,1) 0 dNdx(1,1)  dNdx(3,2) 0 dNdx(1,2) ...
                    dNdx(3,3) 0 dNdx(1,3)  dNdx(3,4) 0 dNdx(1,4) ...
                    dNdx(3,5) 0 dNdx(1,5);
                 0 dNdx(3,1) dNdx(2,1)  0 dNdx(3,2) dNdx(2,2) ...
                    0 dNdx(3,3) dNdx(2,3)  0 dNdx(3,4) dNdx(2,4) ...
                    0 dNdx(3,5) dNdx(2,5)];
        end
    end
    
    %% Pyramid-specific methods
    methods
        function pyramid = fem_pyramid_element(E)
        % Constructor, creates a single or array of elements
        %
        % pyramid = fem_pyramid_element(E)
        %   Creates a set of M pyramid elements given an Mx5 matrix of
        %   node indices, arranged CW w.r.t. the outward normal
        
            if (nargin == 0) 
                pyramid.nodeIdxs = [];
            else
                if (size(E,1)<2)
                    pyramid.nodeIdxs = E;
                else
                    % pre-allocate array
                    pyramid(size(E,1)) = fem_pyramid_element([]);
                    for i=1:size(E,1)
                        pyramid(i).nodeIdxs = E(i,:);
                    end
                end
            end
        end
    end
end

