function t_opt = POCS(t_star, t_centers, R0, R, normA, lambda)
    t = t_star;

    t_min = [-0.5 * normA * lambda, -0.5 * normA * lambda];
    t_max = [ 0.5 * normA * lambda,  0.5 * normA * lambda];
    
    eps_feas = 1e-6;   % feasibility tolerance
    eps_step = 1e-6;   % step tolerance
    K = length(R);
    maxIter = 100;
    
    for iter = 1:maxIter
        t_prev = t;
    
        % Record
    
        %% Projection onto box constraint
        t = min(max(t, t_min), t_max);
    
        %% Projection onto ball constraints
        for k = 1:K
            d = t - t_centers(k,:);
            nd = norm(d,2);
            if nd > R(k)
                t = t_centers(k,:) + R(k) * d / nd;
            end
        end
    
        %% ---------- Feasibility check ----------
        feasible = true;
        
        % (C0) First ball constraint centered at t_star
        if norm(t - t_star,2) > R0 + eps_feas
            feasible = false;
        end
        
        % Box constraint
        if any(t < t_min - eps_feas) || any(t > t_max + eps_feas)
            feasible = false;
        end
        
        % Other ball constraints
        for k = 1:K
            if norm(t - t_centers(k,:),2) > R(k) + eps_feas
                feasible = false;
                break;
            end
        end
        
        % ---------- Stopping criterion ----------
        if feasible && norm(t - t_prev,2) < eps_step
            break;
        end

    end
    t_opt = t;


end