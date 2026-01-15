% Index of the global best particle
[~, gbest] = min(swarm(:,4,1));
gbest_val = swarm(gbest,4,1);
% Updating velocity vectors and positions
for n = 1:n_part
    for d = 1 : dim
        r1 = rand; r2 = rand;
        x = swarm(n,1,d);
        v = swarm(n,2,d);
        p = swarm(n,3,d); %pbest
        g = swarm(gbest,3,d); %gbest
        v = Kap*(v + c1*r1*(p - x) + c2*r2*(g - x));
        v = min(max(v, -vel_clamp(d)), vel_clamp(d));
        x = x + v;
        % Absorbing walls: When a particle hits the boundary of the
        % solution space, the velocity is zeroed in that dimension.
        if x < xmin(d)
            x = xmin(d);
            v = 0; % absorb: kill velocity in this dim
        elseif x > xmax(d)
            x = xmax(d);
            v = 0; % absorb: kill velocity in this dim
        end
        % % Update
        swarm(n,2,d) = v;
        swarm(n,1,d) = x;
    end
end

swarms{i} = swarm;

if i > 1
    b_fitness(i) = min(gbest_val, b_fitness(i - 1));

else
    b_fitness(i) = gbest_val;
end

b_swarm = swarm(gbest,3,:);

display(['Iteration: ' num2str(i) '  Fitness: ' num2str(fitness(gbest)) ' Fitness(best): ' num2str(gbest_val)])