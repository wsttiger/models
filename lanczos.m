function [B,gsv] = lanczos(A,n)

m = size(A,1);
v = rand(m,1);
v = v / norm(v);
vold = zeros(m,1);

% save initial v
vinit = v;

a = zeros(n,1);
b = zeros(n-1,1);
Vsaved = zeros(m,n);
for i = 1:n
    Vsaved(:,i) = v;
    v2 = A*v;
    if (i > 1)
        v2 = v2 - b(i-1)*vold;
    end
    a(i) = v'*v2;
    v2 = v2 - a(i)*v;
    if (i < n)
        b(i) = norm(v2);
        vold = v;
        v = v2 / b(i);
    end
end

B = diag(a) + diag(b,1) + diag(b,-1);

[ev,e] = eig(B);
gsv = Vsaved*ev(:,1);

v = vinit;
gsv = zeros(m,1);
% a = zeros(n,1);
% b = zeros(n-1,1);
for i = 1:n
    % sum for ground state
    gsv = gsv + ev(i,1)*v;
    v2 = A*v;
    if (i > 1)
        v2 = v2 - b(i-1)*vold;
    end
%     a(i) = v'*v2;
    v2 = v2 - a(i)*v;
    if (i < n)
%         b(i) = norm(v2);
        vold = v;
        v = v2 / b(i);
    end
end

gsv = gsv / norm(gsv);


