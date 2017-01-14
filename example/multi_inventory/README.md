An inventory problem example using DP functions. Here the objective is to identify the optimal order for multiple products as a function of the current inventory for each. A key challenge is that there is a space limit in the warehouse. 

Start with the MultiInvDp_scratchpad.m function. 

*Warning* These functions have not been exaustively tested using the latest DP toolbox. So you may encouter some minor errors. More importantly when I developed this, I had access to a whole slew of MATLAB toolboxes, so there may be function calls to functions that were implemented in one of these. In many cases, these functions are pretty simple to write, so we might need to create some of our own.

_Note:_ Because DP exhastively lists all of the states, a key piece of efficient solution techniques is to figure out which aspects of the state space can be safely removed. For example any state that exceeds the capacity of the warehouse can be eliminated. Likewise the set of possible valid decisions for a given starting inventory is limited by remaining warehouse capacity. This state and decision space limiting is handled by the `CombinWithLimits` function. Incedentally, a similar phenomena occurs in the generation expansion planning problem where we know you will never need all of the possible generators of all types, rather there is some planning reserve based upper bound on the the total capacity. The examples for those will take a bit more effort to pull out.

At the time these were written, our (BryanP) research group had not yet figured out how to make ADP work, so there are no comparisons here.

ToDo: Implement with ADP? using as many of the existing functions as possible.
