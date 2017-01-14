function [sample_reward] = inventory_cost(postdec_state, order, demand)

% Ths program calculates current period reward for the current sample.

K = 4;              % per unit fixed cost for placing orders
H = 1;              % per unit holding cost for current inventory

if (order(1) <= 0)
    operating_cost(1) = 0;
else
    operating_cost(1) = K + (2 * order(1));
end

if (order(2) == 0)
    operating_cost(2) = 0;
else
    operating_cost(2) = K + (2 * order(2));
end

operating_cost(3) = K + (2 * order(3));


if (postdec_state(1) < 0)
    postdec_state(1) = 0;
end
    
holding_cost = postdec_state * H;

if (demand > postdec_state)
    demand = postdec_state;
else
end


revenue = 8 * min(demand, postdec_state);

%postdec_state
%demand
%revenue
%holding_cost
%operating_cost

sample_reward  = (revenue - operating_cost - holding_cost);

end
