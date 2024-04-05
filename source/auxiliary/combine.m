function combis = combine(C,n)

expand = 1;
full_expand = C^n;
indices = 1:C;
combis = zeros(full_expand,n);
for k = 1:n
    for kk = 1:length(indices)
        bas = 1+(kk-1)*expand;
        combis(bas:bas+expand-1,k) = indices(kk)*ones(expand,1);
        kkk = k - 1;
        while kkk > 0
            combis(bas:bas+expand-1,kkk) = combis(1:expand,kkk);
            kkk = kkk - 1;
        end
    end
    expand = expand*C;
end