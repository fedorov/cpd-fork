classdef fem_tet_element < fem_element
    %FEM_TET_ELEMENT Linear tetrahedral element
    %
    % Node indices are assumed to be ordered as follows:
    %   1:3, nodes of one face ordered clockwise w.r.t. outward normal
    %   4, final node
    %
    %   Copyright 2013 C. Antonio Sanchez [antonios@ece.ubc.ca]
    
    %% Constants
    properties (Constant, Access=private)
        
        ipntLength = 1;     % number of cubature points
        ipntW = 1.0/6;      % cubature weights
        ipntLoc = [0.5 0.5 0.5];    % cubature point locations
        
        % shape function gradient in natural coordinates
        dN = [-1 1 0 0;
              -1 0 1 0; 
              -1 0 0 1];
          
        % number of nodes in tet elements
        tetNumNodes = 4;
    end
    
    %% Implemented methods
    methods
        
        function n = getNumNodes(~)
            n = fem_tet_element.tetNumNodes;
        end
        
        function N = getN(~,xi,eta,mu)
            N = [1-xi-eta-mu, xi, eta, mu];
        end
        
        function dN = getdN(~,~,~,~)
            dN = fem_tet_element.dN;
        end
        
        function [in] = isNaturalInside(~, xem)
            s = 1-sum(xem,2);
            in = (xem(1) >= 0 && xem(2) >= 0 && xem(3) >= 0 && s >= 0);
        end
        
        function [ipnts, w] = getIPnts(~)
            ipnts = fem_tet_element.ipntLoc;
            w = fem_tet_element.ipntW;
        end
        
    end
    
    %% Over-ridden methods
    methods
        
        function J= getJ(~,~,X) 
            J = fem_tet_element.dN*X;
        end
        
        function dNdx = getShapeGrad(~,~,J)
            dNdx = J\fem_tet_element.dN;
        end
        
        function B = getB(~, dNdx)
            B = [dNdx(1,1) 0 0  dNdx(1,2) 0 0  dNdx(1,3) 0 0 ...
                    dNdx(1,4) 0 0;
                 0 dNdx(2,1) 0  0 dNdx(2,2) 0  0 dNdx(2,3) 0 ...
                    0 dNdx(2,4) 0;
                 0 0 dNdx(3,1)  0 0 dNdx(3,2)  0 0 dNdx(3,3) ...
                    0 0 dNdx(3,4);
                 dNdx(2,1) dNdx(1,1) 0  dNdx(2,2) dNdx(1,2) 0 ...
                    dNdx(2,3) dNdx(1,3) 0  dNdx(2,4) dNdx(1,4) 0;
                 dNdx(3,1) 0 dNdx(1,1)  dNdx(3,2) 0 dNdx(1,2) ...
                    dNdx(3,3) 0 dNdx(1,3)  dNdx(3,4) 0 dNdx(1,4);
                 0 dNdx(3,1) dNdx(2,1)  0 dNdx(3,2) dNdx(2,2) ...
                    0 dNdx(3,3) dNdx(2,3)  0 dNdx(3,4) dNdx(2,4)];
        end
        
        function [ xem ] = getNaturalCoords(~, X, x)
        % Computes the natural coordinates of point x in this element
        %
        % xem = getNaturalCoords(elem, X, x)
        %   Uses a matrix solve to determine the natural coordinates of x
        %
        %   X:  the 4x3 locations of the nodes belonging to this element
        %   x:  the 1x3 spatial coordinates of the point for which to 
        %       determine coordinates
        
            x1 = X(1, :);
            X = X(2:end, :)-repmat(x1,[3,1]);
            x = x-x1;
            
            xem = x/X;
        
        end
        
        function [min, max] = getNaturalCoordBounds(~)
            min = fem_tet_element.naturalLBound;
            max = fem_tet_element.naturalUBound;
        end
        
        function [out, xem] = isInside(this, X, x)
            xem = getNaturalCoords(this, X, x);
            N = getN(this, xem(1), xem(2), xem(3));
            out = (min(N) >= 0 && max(N) <= 1);
        end
        
        function [K, minJ] = getK(this, D, X)
 
            J = getJ(this,[],X); 
            dNdx = getShapeGrad(this,[],J);
            Bj = getB(this,dNdx);
            minJ = det(J);
            K = Bj'*D*Bj*minJ/6;
            
            % more generic:
            % [~, w] = getIPnts(this);
            % n = getNumNodes(this);
            % K = zeros([3*n,3*n]);
            % for i=1:length(w)
            %   J = getJ(this,[],X);
            %   dNdx = getShapeGrad(this,[], J);
            %   B = getB(this,dNdx);
            %   K = K + w(i)*det(J)*B'*D*B;
            % end
        end
    end
    
    %% Tet-specific methods
    methods
        function tet = fem_tet_element(E)
        % Constructor, creates a single or array of elements
        %
        % tets = fem_tet_element(E)
        %   Creates a set of M tetrahedral elements given an Mx4 matrix of
        %   node indices, arranged CW w.r.t. the outward normal
        
            if (nargin == 0) 
                tet.nodeIdxs = [];
            else
                if (size(E,1)<2)
                    tet.nodeIdxs = E;
                else
                    % pre-allocate array
                    tet(size(E,1)) = fem_tet_element([]);
                    for i=1:size(E,1)
                        tet(i).nodeIdxs = E(i,:);
                    end
                end
            end
        end
    end
    
end

