add_bcpd_paths;

%% Hex nodes and elements
  
X_hex_coords = [-1,-1, 0;
                 3,-1, 1;
                 4, 6,-1;
                -2, 5, 0;
                -1, 1, 6;
                 4, 0, 5;
                 3, 7, 5;
                 0, 4, 4];
            
X_hex_coords_multi = [ -1, 1, 0;
                        0, 1, 0;
                        1, 1, 0;
                       -1, 1, 1;
                        0, 1, 1;
                        1, 1, 1;
                       -1, 1, 2;
                        0, 1, 2;
                        1, 1, 2;
                       -1, 0, 0;
                        0, 0, 0;
                        1, 0, 0;
                       -1, 0, 1;
                        0, 0, 1;
                        1, 0, 1;
                       -1, 0, 2;
                        0, 0, 2;
                        1, 0, 2;
                       -1,-1, 0;
                        0,-1, 0;
                        1,-1, 0;
                       -1,-1, 1;
                        0,-1, 1;
                        1,-1, 1;
                       -1, -1, 2;
                        0, -1, 2;
                        1, -1, 2];
            
X_hex_unit = [-1, -1, -1;
               1, -1, -1;
               1,  1, -1;
              -1,  1, -1;
              -1, -1,  1;
               1, -1,  1;
               1,  1,  1;
              -1,  1,  1];
     
E_hex_unit = [ 1 2 3 4 5 6 7 8]; 

E_hex_idxs = [ 1 2 3 4 5 6 7 8];

E_hex_idxs_multi = [ 1  2  5  4 10 11 14 13;
                    10 11 14 13 19 20 23 22;
                     4  5  8  7 13 14 17 16;
                    13 14 17 16 22 23 26 25;
                     2  3  6  5 11 12 15 14;
                    11 12 15 14 20 21 24 23;
                     5  6  9  8 14 15 18 17;
                    14 15 18 17 23 24 27 26 ];  
%% Stiffness parameters

E = 480;
nu = 1/3;

D = fem_material_linear(E,nu);

%% Hex stiffness from objects
elem1 = num2cell(fem_element.create(E_hex_idxs));
elem_multi = num2cell(fem_element.create(E_hex_idxs_multi));

K = fem_stiffness_matrix(X_hex_coords,elem1,D);
K_unit = fem_stiffness_matrix(X_hex_unit, elem1, D);
K_multi = fem_stiffness_matrix(X_hex_coords_multi, elem_multi, D);

%% Comparisons

% load solutions
load fem_hex_stiffness_test_sol;
    
disp('Unit stiffness error:');
disp(max(abs(K_unit(:)-K_hex_unit_sol(:))));

disp('Single hex error:');
disp(max(abs(K(:)-K_hex_single_sol(:))));

disp('Multi hex error:');
disp(max(abs(K_multi(:)-K_hex_multi_sol(:))));

