function A = ILU_mp(A,p)

n = length(A);
LEV = inf(n);
LEV(find(A)) = 0;
LEV(1:n+1:n^2) = 0;    

for i = 2:n   
    r = zeros(1,n);
	for k = 1:i-1   
        if LEV(i,k) > p
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        tmp = A(i,k) * A(k,k+1:n);
        r(k+1:n) = r(k+1:n) + tmp;
        A(i,k+1:n) = A(i,k+1:n) - tmp;
        LEV(i,k+1:n) = min([LEV(i,k+1:n); LEV(i,k) + LEV(k,k+1:n) + 1]);
    end
    zunullen = find(LEV(i,:)>p);
	A(i,zunullen) = 0;
    A(i,i) = A(i,i) - sum(r(zunullen)); 
end
