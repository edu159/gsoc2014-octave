a=load('-ascii','ichol_matrix_10.data');
%a=load('-ascii','ichol_matrix_2000.data');
%a=load('-ascii','ichol_matrix_5000.data');
%a=load('-ascii','matlab_matrix_400.data');
%a=load('-ascii','matlab_matrix_50.data');
%a=load('-ascii','matlab_matrix_10.data');

%a = [3+3i 4-0.2i 4i 3 0; 2-6i 1 0 9+3i 0; 2-2i 0 8+2.3i 9-5.65i 2.2; 0 0 2.3-4.2i 2i 0; 9+2.2i 0 0 2+1i 1-2.34i];

%Este caso da un 0 pivot
%a = [1 4 4 3 0; 2 2 0 9 0; 2 0 8 9 2.2; 0 0 2.3 2 0; 9 0 0 2 1];
%a = [3 4 4 3 0; 2 1 0 9 0; 7 0 8 4 2.2; 0 0 2.3 3 0; 9 0 0 2 1];

%-----------------
%Casos p**
%a = [3 4 4 3 1; 3.2 10 2.3 2 4.5; 9 2 6 2 1;2 0 8 9 2.2; 2 1 9 9 1 ];
%a = [3 4 4 3 1; 3.2 10 2.3 2 4.5; 9 2 6 2 1;2 1 9 9 1; 2 0 8 9 2.2 ];
%** tb

%a = [3 4 4 3 -1.1; 2 0 8 9 2.2; 2 1 9 9 1; 3.2 10 2.3 2 4.5; 9 2 6 2 1];
%a = [3 4 4 3 1; 2 1 9 9 1; 2 0 8 9 2.2; 3.2 9 2.3 2 4.5; 9 2 6 2 1];
%------------------

%a = [3 4 4 3 0; -6 1 0 -3 0; -2 0 8 9 2.2; 0 0 4.2 2 0; 9 0 0 2 2.34];
%a = [3 -4 4 3 1; 2 0 -8 9 1.2; 2 1 9 9 1; 3.2 10 -2.3 2 4.5; 9 2 6 2 1];
%a= [3 -4 4 3 1;  2 1 9 9 1;2 0 -8 9 1.2; 3.2 10 -2.3 2 4.5; 9 2 6 2 1];
%a= [3 1 0 0 4; 3 1 0 0 -2;0 0 8 0 0; 0 4 0 4 -4.5; 0 -1 0 0 1];
% a = [ 13 4 5 9 7 11; 4 11 2 6 0 8; 4 12 2 0 6 0; 8 12 6 16 8 9; 3 12 4 8 12 0; 4  3 8  7  8 9];
%Diag
%a=[0 0 4 3 1; 5 1 2.3 2 4.5; 0 0 0 2 1;0 0 8 0 2.2; 0 0 9 9 1 ];
%ICHOL MATRICES
%a = [10+i 2 3 4 5; 2+2i 20 3 0 2; 2 1 15 0 0; 2 5 8 9 0; 0 0 1 3 9];
%a = [10+2i 2 3 4 5; 0 20 3 0 2; 2 1 15 0 0; 2 5 8 9 0; 0 0 1 3 9];
%a = [10 0 0 0 0; 0 20 0 0 0; 2 1 15 0 0; 2 5 8 9 0; 0 0 1 3 9];
%a = [0.37 -0.05 -0.05 -0.07; -0.05 0.116 0 -0.05; -0.05 0 0.116 -0.05; -0.07 -0.05 -0.05 0.202];
%a = [0.37 -0.04 -0.05 -0.07; -0.04 0.116 0 -0.05; -0.05 0 0.116 -0.05; -0.07 -0.05 -0.05 0.202];
%a = A = sparse ([3 4 4 3 -1.1; 2 0 8 9 2.2; 2 1 9 9 1; 3.2 10 2.3 2 4.5; 9 2 6 2 1]);
a = sparse(a);
%setup.droptol = 0.0
%setup.thresh = 0;
%setup.udiag = 0;
%setup.milu = 'col';
%setup.type = 'ilutp';
setup.type = 'nofill';
setup.michol = 'off';
setup.diagcomp = 0;
setup.droptol = 1e-4;
setup.shape = 'lower';

tic
%L=ICHOL0(a, 'on');
%L=ICHOLT(a, 'on', 0.2);
%s=ilu0(a);
%[L,U] = ILU_crout(a, setup.droptol);

%[L, U, P] = ilu(a, setup);
L = ichol(a, setup);
toc
%L = tril(s,-1) +speye(length(s));
%U = triu(s);
