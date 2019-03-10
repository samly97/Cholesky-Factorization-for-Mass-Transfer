function L = cholesky(A)
  [n,m] = size(A);
  
  if n ~= m
      error('abort! Matrix not square')
  end
  
  if A-1/2*(A+A') ~= zeros(n,n)
      error('matrix not symmetric')
  end
  
  if det(A) < 0
      error('matrix is not positive-definite')
  end
  
  L = zeros(n,n);
  
  for i = 1:n    % column
     for j = i:n % row
         % computing the diagonals
         if j == i
             if j == 1
                 L(j,i) = sqrt(A(j,i));
             else
                 L(j,i) = sqrt(A(j,i) - L(i,1:i-1)*L(i,1:i-1)');
             end
         % elements left of diagonals
         else % j < i
             if i == 1
                L(j,i) = A(j,i)/L(i,i);
             else
                 L(j,i) = (A(i,j)-L(i,1:i-1)*L(j,1:i-1)')/L(i,i);
             end
         end
     end
  end
end