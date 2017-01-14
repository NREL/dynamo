function [sample_reward] = inventory_cost_regress(postdec_state, order, demand)

% Authors:
% Nidhi R. Santen (MIT ESD) and Maria Teresa Pena Alcaraz (IIT, Comillas
% University), with help from Mort Webster (MIT ESD)
% October 2010

% Ths program calculates current period reward for the current sample.

% MAITE: IF WE HAVE K<>0 THEN THE OPTIMAL POLICY IS NOT ORDER TILL OPTIMAL
% STATE
K = 0;              % per unit fixed cost for placing orders
H = 1;              % per unit holding cost for current inventory
UC= 2;              % per unit cost
UR= 8;              % per unit revenue

if (order <= 0)
    operating_cost = 0;
else
    operating_cost = K + (UC * order);
end

% MAITE: This checking should be superfluous
if (postdec_state < 0)
    postdec_state = 0;
end
    
holding_cost = max(postdec_state-demand, 0) * H;
revenue = UR * min(demand, postdec_state);

sample_reward  = (revenue - operating_cost - holding_cost);

end
