function thetas = EM_VmF(z,kappa,mu, n_iter)
%EM_VMF Computes the Expectation Maximization solution for hypersphere 
%   fitting using a von Mises-Fisher prior
%   THETAS = EM_VmF(Z, KAPPA, MU, N_ITER) where Z is a N-by-D matrix of N 
%   vectors of dimension D distributed on a hypersphere, KAPPA>=0 is the 
%   concentration parameter of the von Mise-Fisher prior (KAPPA=0 
%   corresponds to a uniform prior), MU is a vector with norm 1 and is the
%   mean direction of the von Mises-Fisher prior, and N_ITER is the number 
%   of iterations before stopping the algorithm returns a matrix of size 
%   D+1-by-N-ITER where each columns is an iteration of the algorithm. 
%   The first line corresponds to the estimated radius and the next D to 
%   the estimated center and the last to the esitmated noise power.

n=size(z,1);
d = size(z,2);

% Initialization
ct = mean(z);
rt = mean(vecnorm(transpose(z-ct)));
sigma2t = (1/d)*(mean(vecnorm(transpose(z-ct)).^2)-rt^2);
thetas = zeros(d+2,n_iter+1);
thetas(:,1) = [rt;ct';sigma2t];

for k = 1 : n_iter
    
    delta_t = rt*(z-ct) + (sigma2t)*kappa*mu;
    kappa_t = vecnorm(delta_t,2,2)./(sigma2t);
    mu_t = (delta_t)./vecnorm(delta_t,2,2);
    
    alpha_t = besseli(d/2,kappa_t)./besseli(d/2-1,kappa_t).*mu_t;
    alpha_t(isnan(alpha_t(:,1)),:) = mu_t(isnan(alpha_t(:,1)),:); % avoid nans: when bessel function has infinite values, the ratio is taken to 1
    
    f = [sum(sum(alpha_t.*z,2)) ; sum(z)'];
    
    alpha = sum(alpha_t);
    H_inv = 1/(n^2-alpha*alpha')*[[n, -alpha];[-alpha', ((n^2-alpha*alpha')/(n))*eye(d)+(1/n)*(alpha'*alpha)]];
    
    rc = H_inv*f;
    
    rt = rc(1);
    ct = rc(2:end)';
    
    M = sum(vecnorm(z,2,2).^2)-f'*rc;
    
    sigma2t = M/(n*d);
    
    thetas(:,k+1) = [rt;ct';sigma2t];  
    
end

end

