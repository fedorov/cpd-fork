%   CPD_PLOT(X, Y, C); plots 2 data sets. Works only for 2D and 3D data sets.
%
%   Input
%   ------------------ 
%   X           Reference point set matrix NxD;
%   Y           Current postions of GMM centroids;
%   C           (optional) The correspondence vector, such that Y corresponds to X(C,:) 
%
%   See also CPD_REGISTER.

% Copyright (C) 2007 Andriy Myronenko (myron@csee.ogi.edu)
%
%     This file is part of the Coherent Point Drift (CPD) package.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     CPD package is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with CPD package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

function cpd_plot_iter_groupwise(X, Y, T, C)

if nargin<2, error('cpd_plot.m error! Not enough input parameters.'); end;
[m, d]=size(Y);

if d>4, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;
if d<2, error('cpd_plot.m error! Supported dimension for visualizations are only 2D and 3D.'); end;

% for 2D case
if d==2,
%     hold off
%    axis([-1.5 2 -1.5 2]);
 plot(Y(:,1), Y(:,2),'r+');
    hold on
    for i=1:length(X)
        if(i == 1)
            plot(X{i}(:,1), X{i}(:,2),'o', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'b', 'MarkerSize',3)
%             plot(T{i}(:,1), T{i}(:,2),'b.')
        else
            plot(X{i}(:,1), X{i}(:,2),'o', 'MarkerFaceColor','g', 'MarkerEdgeColor', 'g', 'MarkerSize',3)
%             plot(T{i}(:,1), T{i}(:,2),'g.')
        end
    end
    %axis off; 
%     hold off
    axis([-2.5 2.5 -2.5 2.5]);
   
else
    hold off
% for 3D case
   plot3(Y(:,1),Y(:,2),Y(:,3),'r.')
       hold on
    for i=1:length(X)
        if(i == 1)
            plot3(X{i}(:,1), X{i}(:,2), X{i}(:,3),'bo')
        else
            plot3(X{i}(:,1), X{i}(:,2), X{i}(:,3),'go')
        end
    end
    
end
if nargin>3,
    hold on;
    for s=1:length(X)
        if d==2,
%             for i=1:size(Y, 1)
%                 plot([X{s}(C{s}(i),1) Y(i,1)],[X{s}(C{s}(i),2) Y(i,2)]);
%             end
            for i=1:size(T{s}, 1)
                plot([Y(i,1) T{s}(i,1)],[Y(i,2) T{s}(i,2)]);
            end

        else
            for i=1:size(Y, 1)
                plot3([X{s}(C{s}(i),1) Y(i,1)],[X{s}(C{s}(i),2) Y(i,2)],[X{s}(C{s}(i),3) Y(i,3)]);
            end
        end
    end
    hold off;
end

drawnow;