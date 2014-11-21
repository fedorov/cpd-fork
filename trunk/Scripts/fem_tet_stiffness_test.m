add_bcpd_paths;

%% Tet nodes and elements
  
X_tet_coords = [2 3 4;
                6 3 2;
                2 5 1;
                4 3 6];
            
X_tet_coords_multi = [ 0 0 0;
                       1 0 0;
                       1 1 0;
                       0 1 0;
                       0 0 1;
                       1 0 1;
                       1 1 1;
                       0 1 1];
            
X_tet_unit = [0 0 0;
              1 0 0;
              0 1 0;
              0 0 1];
     
E_tet_unit = [ 1 2 3 4 ]; 

E_tet_idxs = [ 1 2 3 4 ];

E_tet_idxs_multi = [ 5 1 4 2;
               7 2 4 3;
               7 6 5 2;
               7 5 8 4;
               7 5 4 2 ];
  
%% Stiffness parameters

E = 480;
nu = 1/3;

D = fem_material_linear.getElasticity(E,nu);
Dmat = fem_material_linear(E, nu);

%% Tet stiffness from objects
elem1 = num2cell(fem_element.create(E_tet_idxs));
elem_multi = num2cell(fem_element.create(E_tet_idxs_multi));

K = fem_stiffness_matrix(X_tet_coords,elem1,D);
K_unit = fem_stiffness_matrix(X_tet_unit, elem1, D);
K_multi = fem_stiffness_matrix(X_tet_coords_multi, elem_multi,D);

%% Tet stiffness method #2
[Phi, V, B] = fem_tet_shape(X_tet_coords, E_tet_idxs);
K2 = sparse(3*size(X_tet_coords,1),3*size(X_tet_coords,1));   

for i=1:size(E_tet_idxs,1)
    Be = B(:,:,i);
    
    Kij = V(i)*Be'*D*Be;
    nodeIdxs = E_tet_idxs(i,:);
    Kidx = 3*(nodeIdxs-1);
        
    Kidx = [Kidx+1;
            Kidx+2;
            Kidx+3];
    K2(Kidx,Kidx) = K2(Kidx,Kidx) + Kij;
end
clear B Be Kidx Kij i;

%% Comparisons

% Unit solution
Ehat = E/(12*(1-2*nu)*(1+nu));
nutilde = 1-2*nu;
nuhat = 1-nu;
K_unit_sol = Ehat*[4-6*nu, 1, 1, -2*nuhat, -nutilde, -nutilde, ...
        -nutilde, -2*nu, 0, -nutilde, 0, -2*nu;
    1, 4-6*nu, 1, -2*nu, -nutilde, 0, -nutilde, -2*nuhat, -nutilde, ...
        0, -nutilde, -2*nu;
    1, 1, 4-6*nu, -2*nu, 0, -nutilde, 0, -2*nu, -nutilde, -nutilde, ...
        -nutilde, -2*nuhat;
    -2*nuhat, -2*nu, -2*nu, 2*nuhat, 0, 0, 0, 2*nu, 0, 0, 0, 2*nu;
    -nutilde, -nutilde, 0, 0, nutilde, 0, nutilde, 0, 0, 0, 0, 0;
    -nutilde, 0, -nutilde, 0, 0, nutilde, 0, 0, 0, nutilde, 0, 0;
    -nutilde, -nutilde, 0, 0, nutilde, 0, nutilde, 0, 0, 0, 0, 0;
    -2*nu, -2*nuhat, -2*nu, 2*nu, 0, 0, 0, 2*nuhat, 0, 0, 0, 2*nu;
    0, -nutilde, -nutilde, 0, 0, 0, 0, 0, nutilde, 0, nutilde, 0;
    -nutilde, 0, -nutilde, 0, 0, nutilde, 0, 0, 0, nutilde, 0, 0;
    0, -nutilde, -nutilde, 0, 0, 0, 0, 0, nutilde, 0, nutilde, 0;
    -2*nu, -2*nu, -2*nuhat, 2*nu, 0, 0, 0, 2*nu, 0, 0, 0, 2*nuhat];

K_unit_sol = sparse(round(K_unit_sol*1e12)/1e12);
clear nuhat nutilde Ehat;

K_sol = [745,  540, 120,  -5,  30,  60,-270, -240,   0,-470, -330,-180;
         540, 1720, 270,-120, 520, 210,-120,-1080, -60,-300,-1160,-420;
         120,  270, 565,   0, 150, 175,   0, -120,-270,-120, -300,-470;
          -5, -120,   0, 145, -90, -60, -90,  120,   0, -50,   90,  60;
          30,  520, 150, -90, 220,  90,  60, -360, -60,   0, -380,-180;
          60,  210, 175, -60,  90, 145,   0, -120, -90,   0, -180,-230;
        -270, -120,   0, -90,  60,   0, 180,    0,   0, 180,   60,   0;
        -240,-1080,-120, 120,-360,-120,   0,  720,   0, 120,  720, 240;
           0,  -60,-270,   0, -60, -90,   0,    0, 180,   0,  120, 180;
        -470, -300,-120, -50,   0,   0, 180,  120,   0, 340,  180, 120;
        -330,-1160,-300,  90,-380,-180,  60,  720, 120, 180,  820, 360;
        -180, -420,-470,  60,-180,-230,   0,  240, 180, 120,  360, 520 ];
    
disp('Unit stiffness error:');
disp(max(abs(K_unit(:)-K_unit_sol(:))));

disp('Single tet error:');
disp(max(abs(K(:)-K_sol(:))));

disp('Estimated error between methods:');
disp(max(abs(K(:)-K2(:))));

load fem_tet_stiffness_test_sol;
disp('Multi tet error:');
disp(max(abs(K_multi(:)-K_tet_multi_sol(:))));