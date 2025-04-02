%{
    ------------------------------------------------------------------------------
    Author(s):    [Vladimir ROKHLIN , Hanwen ZHANG]
    Date:         [April 2024]
    Description:  [Adagpted Gaussian integrator]
    ------------------------------------------------------------------------------
%}

% function [x0,w0] = adapgauss_legpts(N)
%     [x0,w0] = legpts(N);
%     % This routine can be easily replaced with a QR version of O(N^2) work
%     % It is acceptable since we don't need very high order
% end


function [x0,w0] = adapgauss_legpts(N)
    beta = @(n) 0.5*1./sqrt(1-(2*n).^(-2));
    %Jacobi matrix of Legendre
    T = diag(beta(1:(N-1)),-1) + diag(beta(1:(N-1)),+1);
    
    [V,U] = eig(T);
    
    x0 = diag(U);
    w0 = 2*V(1,:).^2;
end
