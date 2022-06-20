del = 0.8:-0.1:0.3;
tau = 0.06:0.02:0.14;
lam = 1.8:0.2:2.0;
rcp = [1 3 4];

count = 0;
for z=0:2
    for i=1:length(del)
        for j=1:length(tau)
            for k=1:length(lam)
                for l=1:length(rcp)
                    g2_rcp_ce_lambda(rcp(l), 1, 2, del(i), tau(j), lam(k), z);
                    %g2_rcp_ce_lambda_halftime(rcp(l), 1, 2, del(i), tau(j), lam(k), z);
                    count = count+1;
                    count/(length(del)*length(tau)*length(lam)*length(rcp)*3)
                end
            end
        end
    end
end