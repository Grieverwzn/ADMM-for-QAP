# ADMM for QAP

Two version of ADMM codes are uploaded. 
- The first one (Source_code (Original QAP model)) is the ADMM algorithm to solve the classic Quadratic Assignment Problem (QAP)
- The second one (Source_code (Household+Office)) is the ADMM algorithm to solve the Quadratic Assignment Problem (QAP) to solving the urban layout problem, where we consider n households and n offices. 

The characteriscs of our method is as follows: 
1. We use a network reformulation to linearize the QAP model. The QAP model is reformulated as a network-flow problem with different types of agents (e.g. builders and householders)
2. The ADMM uses an agent-based block coordination and decompose the primal problem into a series of shortest path subproblems.
3. The augmented Lagrangian method is applied to generated feasible solutions efficiently. We also provide a method to linearize the quadratic augmented term in the mixed integer programmign model.

