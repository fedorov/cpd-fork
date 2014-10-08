add_bcpd_paths;

%% Load prostate and make another instance

load SSM;
m = SSM.mean;

n = m + reshape(SSM.mods(:,1:4)*[2;1;1;1],[size(m,2),size(m,1)])';
n = n(1:3:end,:);

n = n+repmat([4,4,4],[size(n,1),1]);

%% Perform CPD in a loop, animating

tn = n;
sigma2=0.2;

clf;
plot3(m(:,1), m(:,2), m(:,3),'.r','MarkerSize',10);
axis([-5,5,-5,5,-5,5]);
hold on;
plot3(tn(:,1), tn(:,2), tn(:,3),'ob','MarkerSize', 3);
axis([-5,5,-5,5,-5,5]);
pause(0.1);

for i=1:50
    
    [tn, R, t, s, P, sigma2] = cpd_rigid(m, tn, 0, 1e-10, 1, [],[],[], sigma2);
    
    clf;
    plot3(m(:,1), m(:,2), m(:,3),'.r','MarkerSize',10);
    axis([-5,5,-5,5,-5,5]);
    hold on;
    plot3(tn(:,1), tn(:,2), tn(:,3),'ob','MarkerSize', 3);
    axis([-5,5,-5,5,-5,5]);
    pause(0.05);
end