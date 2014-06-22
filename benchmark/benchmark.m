a=load('-ascii','matlab_matrix_5000.data');

%a=load('-ascii','matlab_matrix_2000.data');
printf("acabo");
fflush(stdout);
%a = a(1:9,1:9);
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

%a = [ 13 4 5 6 7 6; 2 11 0 6 0 8; 4 12 2 0 6 0; 8 4 6 16 8 9; 3 2 4 8 12 0; 4  1 1  7  8 9] 
%Diag
%a=[0 0 4 3 1; 5 1 2.3 2 4.5; 0 0 0 2 1;0 0 8 0 2.2; 0 0 9 9 1 ];
a = sparse(a);
setup.droptol = 0.2;
setup.thresh = 0;
setup.udiag = 0;
setup.milu = 'off';

setup.type = 'crout';
tic
%s=ilu0(a);
%[L,U] = ILU_crout(a, setup.droptol);
%[s, P]=ILU_pc(a,setup.droptol, setup.thresh );
%s = sparse(s);
[L, U, P] = ilu(a, setup);
toc
%L = tril(s,-1) +speye(length(s));
%U = triu(s);
