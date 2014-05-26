%a=load('-ascii','matlab_matrix_50.data');
%a=load('-ascii','matlab_matrix_10.data');
%a = [3+3i 4-0.2i 4i 3 0; 2-6i 0-0i 0 9+3i 0; 2-2i 0 8+2.3i 9-5.65i 2.2; 0 0 2.3-4.2i 2i 0; 9+2.2i 0 0 2+1i 1-2.34i];
%a = [3 4 4 3 0; 2 1 0 9 0; 2 0 8 9 2.2; 0 0 2.3 2 0; 9 0 0 2 1];
%a = [3 4 4 3 0; 2 1 0 9 0; 7 0 8 4 2.2; 0 0 2.3 3 0; 9 0 0 2 1];
a = [3 4 4 3 1; 2 1 9 9 1; 2 0 8 9 2.2; 3.2 10 2.3 2 4.5; 9 2 6 2 1];
a = sparse(a);
tic
%s=ilu0(a);
%s=ILU_0(a);
droptol = 0.1;
thresh = 0.5
[s, P]=ILU_pc(a,droptol, thresh );
toc
s = sparse(s);
L = tril(s,-1) +speye(length(a));
U = triu(s);
