%a=load('-ascii',);
a = [3+3i 4-0.2i 4i 3 0; 2-6i 0-0i 0 9+3i 0; 2-2i 0 8+2.3i 9-5.65i 2.2; 0 0 2.3-4.2i 2i 0; 9+2.2i 0 0 2+1i 1-2.34i];
a = sparse(a);
tic
s=ilu0(a);
%s=ILU_0(a);
toc
L = tril(s,-1) +speye(length(a));
U = triu(s);
