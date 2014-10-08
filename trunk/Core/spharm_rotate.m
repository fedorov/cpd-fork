function rvec = spharm_rotate(alpha, beta, gamma, fvec, max_degree)
sgm = {};

for d=1:max_degree
    n = 2*d+1;
    M = (-1).^(1:n^2);
    sgm{d} = reshape(M,n,n);
end

%beta  = -beta; % caused by nTHETA = pi/2-THETA (original nTHETA = pi/2+THETA)?

fnum = (max_degree+1)^2;
for k=1:length(beta)
    D{k} = sparse(fnum,fnum);
end
for l = 0:max_degree
    M = [];
    for m = -l:l     % actual order (-l => l)
        ma = m+l; % l^2+ma+1
        for n = -m:l % actual order (-l => l) 
            d_beta = zeros(size(beta));
            na = n+l; % l^2+na+1
			for t = max(0,n-m):min(l+n,l-m)
              d_beta_temp = (-1)^(t)*(sqrt(factorial(l+n)*factorial(l-n)*factorial(l+m)*factorial(l-m))) ...
                            /(factorial(l+n-t)*factorial(l-m-t)*factorial(t+m-n)*factorial(t));
              d_beta_temp = d_beta_temp * (cos(beta/2)).^(2*l+n-m-2*t) .* (sin(beta/2)).^(2*t+m-n);
              d_beta = d_beta + d_beta_temp;
			end
            % use symmetry: d(l,m,n,beta) = d(l,-n,-m,beta); d(l,m,n,beta) = (-1)^(m+n)d(l,n,m,beta)
            M(ma+1,na+1,:) = exp(-i*m*alpha).*d_beta.*exp(-i*n*gamma); 
        end
    end 
    idx = (l^2+1):(l+1)^2; 
    for k=1:length(beta)
        if l>0
            Mx = M(:,end:-1:1,k).*(1-eye(size(M(:,:,k)))); Mx = Mx(end:-1:1,:);
            M(:,:,k) = M(:,:,k) - real(Mx).*sgm{l} + i*imag(Mx).*sgm{l};
        end       
        D{k}(idx,idx) = M(:,:,k); 
    end
end
for k=1:length(beta)
    rvec(:,:,k) = D{k}*fvec;
end

return;
