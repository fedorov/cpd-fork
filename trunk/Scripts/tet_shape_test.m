add_bcpd_paths;

%% Nodes and elements

X = [ 0 0 0;
      1 0 0;
      1 1 0;
      0 1 0;
      0 0 1;
      1 0 1;
      1 1 1;
      0 1 1 ];

E = [ 5 1 4 2;
      7 2 4 3;
      7 6 5 2;
      7 5 8 4;
      7 5 4 2 ];
  

%% Compute shape functions

[Phi, V, B] = fem_tet_shape(X, E);

disp('Total volume: ');
disp(sum(V));

disp('Evaluating shape functions at nodes...');
for i=1:size(E,1)
   disp(['Element ',num2str(i)]);
   N = X(E(i,:),:);
   mid = sum(N,1)/4;
   % get all nodes plus centre
   A = [ones(size(E,2),1) N; [1 mid]];    
   P = Phi(:,:,i);          % get all shape functions
   disp('Expecting: [I; 0.25 0.25 0.25 0.25]');
   disp(A*P');
end