function X = CGLS(A,b,NK,reorth)
%CGLS Conjugate gradient algorithm applied implicitly to the normal equations.
%
% [X,rho,eta,F] = cgls(A,b,k,reorth,s)
%
% Performs k steps of the conjugate gradient algorithm applied
% implicitly to the normal equations A'*A*x = A'*b.
%
% The routine returns all k solutions, stored as columns of
% the matrix X.  The corresponding solution and residual norms
% are returned in the vectors eta and rho, respectively.
%
% If the singular values s are also provided, cgls computes the
% filter factors associated with each step and stores them
% columnwise in the matrix F.
%
% Reorthogonalization of the normal equation residual vectors
% A'*(A*X(:,i)-b) is controlled by means of reorth:
%    reorth = 0 : no reorthogonalization (default),
%    reorth = 1 : reorthogonalization by means of MGS.

% References: A. Bjorck, "Numerical Methods for Least Squares Problems",
% SIAM, Philadelphia, 1996.
% C. R. Vogel, "Solving ill-conditioned linear systems using the
% conjugate gradient method", Report, Dept. of Mathematical
% Sciences, Montana State University, 1987.

% Per Christian Hansen, IMM, July 23, 2007.

% Modified by Benjamin Shearer (2025).


[~,n] = size(A);
% X = zeros(n,NK);

% Initialization.
if (reorth==1)
    ATr = zeros(n,NK+1);
end

% Prepare for CG iteration.
x = zeros(n,1);   % No guess
% x = b;              % Initialize as b  
d = A'*b;
r = b;
normr2 = d'*d;
if (reorth==1)
    ATr(:,1) = d/norm(d);
end

% % Make coordinate system:
% yy = P.ds_sam*(-(P.N/2) : (P.N/2-1));
% xx = yy';
% rr = sqrt(xx.^2 + yy.^2);
% 
% % Total power:
% pwr = sum(b.^2);

% figure;

% Calculate error of assuming m=b:
error = zeros(1,NK+1);
error(1) = norm(A*b - b);

% Initialize solution matrix:
xsol = zeros(n,NK+1);
xsol(:,1) = b;

% Iterate.
nk = 1;
solutionFound = 0;
while nk <= NK && solutionFound == 0

    % Update x and r vectors.
    Ad = A*d;
    alpha = normr2/(Ad'*Ad);
    x = x + alpha*d;

    % Calculate error for this iteration:
    error(nk+1) = norm(A*x - b);

    % If error is greater than prior iteration, set sol as previous iter:
    if (error(nk+1) > error(nk)) && nk > 1
        X = xsol(:,nk);
        solutionFound = 1;
%         fprintf('kopt = %d. %.2e to %.2e\n',nk-1,error(nk+1),error(nk))
        fprintf('kopt = %d\n',nk-1)

    end

    x = abs(x);
    xsol(:,nk+1) = x;
    
    r  = r - alpha*Ad;
    s  = A'*r;
    
    % Reorthogonalize s to previous s-vectors, if required:
    if (reorth==1)
        for i=1:nk
            s = s - (ATr(:,i)'*s)*ATr(:,i);
        end
        ATr(:,nk+1) = s/norm(s);
    end
    
    % Update d:
    normr2_new = s'*s;
    beta = normr2_new/normr2;
    normr2 = normr2_new;
    d = s + beta*d;

    % Increment step:
    nk = nk + 1;
end

% If last iteration has lowest error, set that as the solution:
if solutionFound == 0
    X = xsol(:,end);
%     fprintf('kopt = %d\n',nk-1)
end