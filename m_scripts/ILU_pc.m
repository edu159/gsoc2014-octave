function A = ILU_pc(A, tau)

B = A;
n = length(A);

for i = 2:n
	for k = 1:i-1
        if ( (abs(A(i,k)) < tau) && (B(i,k)==0) )
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        A(i,k+1:n) = A(i,k+1:n) - A(i,k) * A(k,k+1:n);
    end
	A(i,find( (abs(A(i,:)) < tau) & (B(i,:)==0) )) = 0;
end

 
