function A = ILU_m(A)

B = A;
n = length(A);

for i = 2:n   
	r = 0;
	for k = 1:i-1   
        if B(i,k)==0
            continue
        end
        A(i,k) = A(i,k) / A(k,k);
        for j = k+1:n
            if B(i,j)==0
                r = r + A(i,k) * A(k,j);
		printf("tl: %f ", A(i,k));
		printf("data: %f \n", A(k,j));
                continue
            end
            A(i,j) = A(i,j) - A(i,k) * A(k,j); 
        end
    end

    A(i,i) = A(i,i) - r;    
end
