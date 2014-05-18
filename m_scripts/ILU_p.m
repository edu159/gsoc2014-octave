function A = ILU_p(A,p)

n = length(A);
LEV = inf(n);
LEV(find(A)) = 0;
LEV(1:n+1:n^2) = 0;  

for i = 2:n   
	for k = 1:i-1  
        if LEV(i,k) > p
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
        LEV(i,k+1:n) = min([LEV(i,k+1:n); LEV(i,k) + LEV(k,k+1:n) + 1]);
    end
	A(i,find(LEV(i,:)>p)) = 0; 
end
