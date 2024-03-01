# ClustEnergy OpTool
ClustEnergy OpTool version 1.0.0 is a useful tool for evaluating the energy flexibility potential offered by an aggregate of buildings (e.g., cluster of buildings), which is useful for resource planning in future scenarios. 
It is a simple Python tool based on linear programming that allows simulation of different demand management strategies in different user-defined clusters of buildings. 
Through archetype-based approach, it is possible to define a representative cluster of buildings subjected to a demand management strategy during a simulation period and reference location. 
Then, from the comparison with a baseline scenario (BL), it is possible to study the flexible behavior of a cluster of buildings subjected to a demand response (DR) strategy. 
Specifically, through centralized optimal control, a peak shaving strategy or a load shifting strategy can be implemented. 
Although the currently available archetypes refer to the Italian case study, the tool given its simplicity is easily modified and adaptable to the needs of the individual use.  

For further information on ClustEnergy OpTool operation, refer to the documentation.
For more information on the methodology, see references below.

# Dependencies
ClustEnergy OpTool is implemented in Python 3.9 and uses the following libraries:
- pandas
- numpy
- scipy
- pvlib
- matplotlib
- os
- sys
- psychrolib
- math
- datetime

# Contacts
- Patricia Ercoli p.ercoli@pm.univpm.it
- Alice Mugnini a.mugnini@univpm.it
- Alessia Arteconi a.arteconi@univpm.it

# Revision history
29 February 2024, first version

# References
- Mugnini, A., Polonara, F., Arteconi, A., 2021. Energy Flexibility of Clusters of Buildings: development of an assessment tool. Proceedings of Building Simulation 2021: 17th Conference of IBPSA, IBPSA, 167-174.
- Ercoli, P., Mugnini, A., Caresana, F., Arteconi, A., 2023. Flexible heat pumps in clusters of buildings: energy flexibility quantification of space cooling loads. Proceedings of the 26th IIR International Congress of Refrigeration: Paris, IIR.
- Ercoli, P., Mugnini, A., Polonara, F., Arteconi, A., 2023. Flexible cooling demand in cluster of buildings: energy flexibility quantification in presence of central or distributed photovoltaic generation. Proceedings of Building Simulation 2023: 18th Conference of IBPSA, IBPSA, 3453-3460.
