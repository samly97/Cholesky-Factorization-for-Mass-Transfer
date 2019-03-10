function exploratoryprob
    P = 0.1; % atm
    R = 16/195; % atm*L/mol/K
    T = 27 + 273.15; % K
    Cas = P/(R*T)/1000; % mol/cm^3
    
    hi = 4/25;
    hj = 6/25;
    xn = 0 + hi:hi:4 - hi; % internal nodes
    n = length(xn); % # of int nodes in across row/column

    A = zeros(n^2,n^2);
    b = zeros(n^2,1);
    
    for i = 1:n % x
        for j = 1:n % y
            r = n*(j-1)+i; % row/eqn
            c = r; % column
            
            % dealing with corner nodes
            if (i == 1 || i == n) && (j == 1 || j == n)
                % (1,1) or (n,n)
                if i == j
                    % (1,1)
                    if i == 1
                        A(r,c) = -2*(hi^2+hj^2);
                        A(r,c+1) = hj^2;
                        A(r,c+n) = hi^2;
                    % (n,n)
                    else
                        A(r,c-n) = hi^2;
                        A(r,c-1) = hj^2;
                        A(r,c) = -2*(hi^2+hj^2);
                        b(r,1) = -hi^2*Cas;
                    end
                % (1,n) or (n,1)
                else
                    % (n,1)
                    if i == n
                        A(r,c-1) = hj^2;
                        A(r,c) = -2*(hi^2+hj^2);
                        A(r,c+n) = hi^2;
                    % (1,n)
                    else
                        A(r,c-n) = hi^2;
                        A(r,c) = -2*(hi^2+hj^2);
                        A(r,c+1) = hj^2;
                        b(r,1) = -hi^2*Cas;
                    end
                end
            % nodes (1, 2-(n-1))
            elseif i == 1
                A(r,c-n) = hi^2;
                A(r,c) = -2*(hi^2+hj^2);
                A(r,c+1) = hj^2;
                A(r,c+n) = hi^2;
            % nodes (n, 2-(n-1))
            elseif i == n
                A(r,c-n) = hi^2; 
                A(r,c-1) = hj^2;
                A(r,c) = -2*(hi^2+hj^2);
                A(r,c+n) = hi^2;
            % nodes (2-(n-1), 1)
            elseif j == 1
                A(r,c-1) = hj^2;
                A(r,c) = -2*(hi^2+hj^2);
                A(r,c+1) = hj^2;
                A(r,c+n) = hi^2;
            % nodes (2-(n-1), n)    
            elseif j == n
                A(r,c-n) = hi^2;
                A(r,c-1) = hj^2;
                A(r,c) = -2*(hi^2+hj^2);
                A(r,c+1) = hj^2;
                b(r,1) = -hi^2*Cas;
            % every other interior node
            else
                A(r,c-n) = hi^2;
                A(r,c-1) = hj^2;
                A(r,c) = -2*(hi^2+hj^2);
                A(r,c+1) = hj^2;
                A(r,c+n) = hi^2;
            end
        end
    end
    
    A=-A; b=-b;
    
    t = cputime;
    L = cholesky(A);
    e = cputime-t;
    d = L\b;
    ca = L'\d*100^3; % mol/m^3
    
    Ca = zeros(length(0:hi:4),length(0:hj:6));
    Ca(2:length(0:hi:4)-1,length(0:hi:4)) = Cas*100^3;
    
    for i = 1:n
        Ca(2:length(0:hi:4)-1,1+i) = ca((i-1)*n+1:i*n,1);
    end
    
    [c,h] = contour((0:hi:4)', (0:hj:6)', Ca');
    clabel(c,h)
    xlabel('x (cm)'), ylabel('y (cm)')
    tit = sprintf('Contour Plot (Cholesky %d Nodes - %.2f Seconds).png', n, e);
    saveas(gcf,tit)
    
    e1 = cputime-e;
    [L, U] = lu(A);
    e2 = cputime-e1;
    
    Ca = zeros(length(0:hi:4),length(0:hj:6));
    Ca(2:length(0:hi:4)-1,length(0:hi:4)) = Cas*100^3;
    d = L\b;
    ca_lu = U\d*100^3;
    
    for i = 1:n
        Ca(2:length(0:hi:4)-1,1+i) = ca_lu((i-1)*n+1:i*n,1);
    end
    
    [c,h] = contour((0:hi:4)', (0:hj:6)', Ca');
    clabel(c,h)
    xlabel('x (cm)'), ylabel('y (cm)')
    tit = sprintf('Contour Plot (LU %d Nodes - %.2f Seconds).png', n, e2);
    saveas(gcf,tit)
end