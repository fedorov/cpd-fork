clear,close all hidden, 
disp('exemple1 ....');


d          = 3;
Nx         = 10;
Ny         = 100;
x          = randn(d , Nx);
y          = randn(d , Ny);
w          = rand(1 , Nx);
h          = 2;

e          = 10;
p          = 6;
K          = 5;

v1         = dval(x , y , w , h);

[xc , A_k] = fgt_model(x , w , h , e , K , p);

v2         = fgt_predict(y , xc , A_k , h , e);


disp(sprintf('error = %5.4f',norm(v1 - v2)));

pause,
clear,close all hidden, 
disp('exemple2 ....');



d          = 2;
R          = [2 , 0.4 ; 0.4  3];
Nx         = 1000;

h          = 1;
e          = 10;
K          = round(sqrt(Nx));
p          = 6;



vect       = (-5:0.3:5);
Ny         = length(vect);
w          = (1/Nx)*ones(1 , Nx);

x          = (chol(R)'*randn(d , Nx));

[X , Y]    = meshgrid(vect);
y          = [X(:) , Y(:)]';

[xc , A_k] = fgt_model(x , w , h , e , K , p);
vy         = fgt_predict(y , xc , A_k , h , e , K , p);

densite    = reshape( vy , Ny , Ny);

figure
set(gcf , 'renderer' , 'opengl');
surfc(X , Y , densite)
shading interp
lighting phong

light
alpha(0.5);
hold on
h = plot(x(1 , :) , x(2 , :) , 'r+' , xc(1 , :) , xc(2 , :) , 'ko' , 'markersize' , 10);
hold off
legend(h(1:2) , '\bf{x}' , '\bf{x}_{fgt}' , 0);
view(2)
colorbar
