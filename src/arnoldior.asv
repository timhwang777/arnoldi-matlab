load west0479
A = west0479;
%k = 20;

% exact eigenvalue
lam = eig(full(A));
figure(2)
hold on
plot(real(lam),imag(lam),'r+');

idx = 1;
%% start calculate range k
for k = 10:10:100
n = length(A);
V = zeros(n,k); % orthonormal basis
H = zeros(k,k); % upper hessenberg matrix
v = ones(n,1);

V(:,1) = v/norm(v);

% Gram-Schmidt
for j = 1:k
V(:,j+1) = A*V(:,j); % compute w

for i = 1:j
H(i,j) = V(:,i)'*V(:,j+1);
V(:,j+1) = V(:,j+1)-H(i,j)*V(:,i);
end
% normalization
H(j+1,j) = norm(V(:,j+1));

% -->reorthogonalization process HERE?
% without any if statement
for l = 1:j
mu = V(:,l)'*V(:,j+1);
V(:,j+1) = V(:,j+1)-V(:,l)*mu;
H(l,j) = H(l,j) + mu;
end
H(j+1, j) = norm(V(:,j+1));

V(:,j+1) = V(:,j+1)/H(j+1,j);
% compute residual
Ra(j,idx) = norm(A*V(:,j) - V(:,j+1)*H(:,j)');
Rb(j,idx) = norm(speye(k) - V(:,j+1)'*V(:,j+1));
end

H(k+1,:) = []; 
Rz = eig(full(H));
plot(real(Rz),imag(Rz),'bo');

% calculate relative errors
lumk = max(lam);
muk = max(Rz);
relerror(idx,:) = abs(lumk - muk) / abs(lumk);
idx = idx + 1;
%% end of the range k
end
hold off