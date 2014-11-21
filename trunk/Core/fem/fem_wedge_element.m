classdef fem_wedge_element < fem_element
    %FEM_WEDGE_ELEMENT Linear wedge element
    %
    % Node indices are assumed to be ordered as follows:
    %   1:3, nodes of one triangular face ordered clockwise 
    %        w.r.t. outward normal
    %   4:6, corresponding nodes on opposite face (will be CCW)
    % 
    %   Copyright 2013 C. Antonio Sanchez [antonios@ece.ubc.ca]
    
    %% Constants
    properties (Constant, Access=private)
        
        INTEGRATION_COORDS_GAUSS_6 = ...
           [ 1/6, 1/6, -1/sqrt(3), 1/6;
             4/6, 1/6, -1/sqrt(3), 1/6;
             1/6, 4/6, -1/sqrt(3), 1/6;
             1/6, 1/6,  1/sqrt(3), 1/6;
             4/6, 1/6,  1/sqrt(3), 1/6;
             1/6, 4/6,  1/sqrt(3), 1/6];
        
        
        % cubature point locations
        ipntLength = 6;     % number of cubature points
        % cubature weights
        ipntW = fem_wedge_element.INTEGRATION_COORDS_GAUSS_6(:,4)';  
        % cubature point locations
        ipntLoc = fem_wedge_element.INTEGRATION_COORDS_GAUSS_6(:,1:3);    
        
        % number of nodes in wedge elements
        wedgeNumNodes = 6;
    end
    
    %% Implemented methods
    methods
        function n = getNumNodes(~)
            n = fem_wedge_element.wedgeNumNodes;
        end
        
        function N = getN(~,xi,eta,mu)
            mm = (1-mu);
            mp = (1+mu);
            xes = (1-xi-eta);
            
            N = [xes.*mm, xi.*mm, eta.*mm, xes.*mp, xi.*mp, eta.*mp]/2;
        end
        
        function dN = getdN(~,xi,eta,mu)
            mm = (1-mu);
            mp = (1+mu);
            xes = (1-xi-eta);
            
            dN = [ -mm,  mm,    0, -mp,  mp,   0;
                   -mm,   0,   mm, -mp,  0,   mp;
                  -xes, -xi, -eta, xes,  xi, eta]/2;
        end
        
        function [in] = isNaturalInside(~, xem)
           s = xem(1)+xem(2);
           in = (xem(1) >= 0 && xem(1) <= 1 && xem(2) >= 0 && ...
                 xem(2) <= 1 && xem(3) >= -1 && xem(3) <= 1 && ...
                 s >= 0 && s <= 1 );
        end
        
        function [ipnts, w] = getIPnts(~)
            ipnts = fem_wedge_element.ipntLoc;
            w = fem_wedge_element.ipntW;
        end
        
    end
    
    %% Over-ridden methods
    methods
         function B = getB(~, dNdx)
            B = [dNdx(1,1) 0 0  dNdx(1,2) 0 0  dNdx(1,3) 0 0 ...
                    dNdx(1,4) 0 0  dNdx(1,5) 0 0  dNdx(1,6) 0 0;
                 0 dNdx(2,1) 0  0 dNdx(2,2) 0  0 dNdx(2,3) 0 ...
                    0 dNdx(2,4) 0  0 dNdx(2,5) 0  0 dNdx(2,6) 0;
                 0 0 dNdx(3,1)  0 0 dNdx(3,2)  0 0 dNdx(3,3) ...
                    0 0 dNdx(3,4)  0 0 dNdx(3,5)  0 0 dNdx(3,6);
                 dNdx(2,1) dNdx(1,1) 0  dNdx(2,2) dNdx(1,2) 0 ...
                    dNdx(2,3) dNdx(1,3) 0  dNdx(2,4) dNdx(1,4) 0 ...
                    dNdx(2,5) dNdx(1,5) 0  dNdx(2,6) dNdx(1,6) 0;
                 dNdx(3,1) 0 dNdx(1,1)  dNdx(3,2) 0 dNdx(1,2) ...
                    dNdx(3,3) 0 dNdx(1,3)  dNdx(3,4) 0 dNdx(1,4) ...
                    dNdx(3,5) 0 dNdx(1,5)  dNdx(3,6) 0 dNdx(1,6);
                 0 dNdx(3,1) dNdx(2,1)  0 dNdx(3,2) dNdx(2,2) ...
                    0 dNdx(3,3) dNdx(2,3)  0 dNdx(3,4) dNdx(2,4) ...
                    0 dNdx(3,5) dNdx(2,5)  0 dNdx(3,6) dNdx(2,6)];
        end
    end
    
    %% Wedge-specific methods
    methods
        function wedge = fem_wedge_element(E)
        % Constructor, creates a single or array of elements
        %
        % wedge = fem_wedge_element(E)
        %   Creates a set of M wedge elements given an Mx6 matrix of
        %   node indices, arranged CW w.r.t. the outward normal
        
            if (nargin == 0) 
                wedge.nodeIdxs = [];
            else
                if (size(E,1)<2)
                    wedge.nodeIdxs = E;
                else
                    % pre-allocate array
                    wedge(size(E,1)) = fem_wedge_element([]);
                    for i=1:size(E,1)
                        wedge(i).nodeIdxs = E(i,:);
                    end
                end
            end
        end
    end
end

