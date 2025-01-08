function xsol1 = BiCGSTAB(A,b,NK)
%BICGSTAB Stabilized biconjugate gradient method.
%==========================================================================

% Initialize:
N2 = length(b);
x = zeros(size(b));
r = b - A*x;
% R = r;
% R = ones(size(r));
R = flip(r);
% R = rand(size(b));
p = r;

% Calculate error of assuming m=b:
error = zeros(1,NK+1);
error(1) = norm(A*b - b);

% Initialize solution matrix:
xsol = zeros(N2,NK+1);
xsol(:,1) = b;

% figure;

% Iterate...
solutionFound = 0;
nk = 1;
while nk <= NK && solutionFound == 0

%     imagesc(reshape(sqrt(xsol(:,nk)),[sqrt(N2),sqrt(N2)])),colorbar;axis square
%     title(num2str(nk))
%     drawnow
%     pause(1)

    Ap = A*p;
    alpha = dot(r,R)/dot(Ap,R);
    s = r - alpha*Ap;
    As = A*s;
    omega = dot(As,s)/dot(As,As);
    x = x + alpha*p + omega*s;
%     x = abs(x);
    r_old = r;
    r = s - omega*As;
    beta = (alpha/omega)*dot(r,R)/dot(r_old,R);
    p = r + beta*(p - omega*Ap);

    % Calculate error for this iteration:
    error(nk+1) = norm(A*x - b);

    % If error is greater than prior iteration, set sol as previous iter:
    if (error(nk+1) > error(nk)) && nk > 1
        xsol1 = xsol(:,nk);
        solutionFound = 1;
        fprintf('kopt = %d\n',nk-1)
    end

    x = abs(x);
    xsol(:,nk+1) = x;

    nk = nk + 1;
end

% If last iteration has lowest error, set that as the solution:
if solutionFound == 0
    xsol1 = xsol(:,end);
    fprintf('kopt = %d\n',nk-1)
end

end

