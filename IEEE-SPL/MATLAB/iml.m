function thetas = iml(z, n_iter)
%IML Computes the Iterative Maximum Likelihood solution for hypersphere 
%   fitting
%   THETAS = iml(Z, N_ITER) where Z is a N-by-D matrix of N vectors of
%   dimension D distributed on a hypersphere and N_ITER is the number of 
%   iterations before stopping the algorithm returns a matrix of size D+1-
%   by-N-ITER where each columns is an iteration of the algorithm. 
%   The first line corresponds to the estimated radius and the others to 
%   the estimated center.
%
%   W. Li, J. Zhong, T. A. Gulliver, B. Rong, R. Hu, and Y. Qian, “Fitting 
%   Noisy Data to a Circle: A Simple Iterative Maximum Likelihood Approach,” 
%   in Proc. IEEE Int. Conf. on Communications, July 2011, pp. 1 – 5. 
n=size(z,1);
d = size(z,2);

ct = mean(z);

thetas = zeros(d+1,n_iter);

for k = 1 : n_iter
    rt = mean(vecnorm(z-ct,2,2));
    ct = mean(z+rt*(ct-z)./vecnorm(ct-z,2,2));
    thetas(:,k) = [rt;ct'];
end

end

