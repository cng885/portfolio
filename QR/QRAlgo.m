function [lambda, vkplus1] = QRAlgo( a )
% This function implements the QR method
%  to solve for the eigenvalues of a square 
% matrix, a. For a given error tolerance, when
% we form an upper triangular matrix that contains
% our eigenvalues on the diagonal.
% 
%@param a a square nxn matrix
%@return lambda the eigenvalues
%@return vkplus1 the corresponding eigenvectors
%@author Chase Ginther
%@date 2016.09.27
    ak = a;
    error = 1;
    k = 0;
    vk = eye(size(a));
    
    while error > 1e-9
        
        [qk, rk] = qr(ak);
        
        akplus1 = rk*qk;
        
                
        vkplus1 = vk * qk;
        
        error = norm(diag(akplus1)-diag(ak));
        
        ak = akplus1;
        vk = vkplus1;
        k = k + 1;
    end
    
    lambda = diag(ak);

end

