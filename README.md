# fastMethod

This repository contains the code and computed gMCSs referenced in the paper 'A fast method for extracting essential and synthetic
lethality genes in GEM models'

The file FM.py contains the code of the algorithm and is located in the 'code' directory.  This module requires:
- The Python modules cobra, numpy, more_itertools, multiprocessing and gc
- A running cplex library or gurobi under Python must be installed in order to run the notebooks.

In order to compute all gMCSs up to a given length for a model, the following steps must be followed:
- Import the model
- Import the module FM by writing 'import FM'
- Set the variable maxLength as the maximum length of the desired gMCSs
- Execute the order 'gMCSs=createGMCS(_model,_maxLength,_solver="gurobi")'

The necessary parameters, _model and _maxLength, are the model and maxLength defined. Finally the parameter _cplex allows the user to changed the solver to be used.

In the 'Notebooks' directory there are two notebooks including examples on the use of the module to compute gMCSs and genetic interventions to ensure the coupled growth of a desired bioproducts.

We have also included in the directories 'models' and 'gMCSs' the models used in the paper and their corresponding gMCSs. The model for Human-GEM and the sets of gMCSs are compressed using Pickle.



