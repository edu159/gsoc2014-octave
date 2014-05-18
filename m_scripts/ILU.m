function A = ILU(A)

B = A;
n = length(A);

for k = 1:n-1   
	for i = k+1:n   
        if B(i,k)==0
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        for j = k+1:n
            if B(i,j)==0
                continue
            end
            A(i,j) = A(i,j) - A(i,k) * A(k,j); 
        end
    end
end
