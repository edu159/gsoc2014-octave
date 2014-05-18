function A = ILU_0(A)

B = A;
n = length(A);

for i = 2:n   
	for k = 1:i-1   
        if B(i,k)==0
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        for j = k+1:n
            if B(i,j)==0
                continue
            end
            A(i,j) = A(i,j) - A(i,k) * A(k,j);
	    printf("i,k,j : %d,%d,%d", i, k,j)
	    A
        end

    end
end
