This folder contains 4 files apart from README.txt, and a sub-folder. 

1. 'Neant_Kissler_Wolfel_input_monolix.csv' is the digitized viral dynamics data from three papers that is input into monolix.
2. 'vdm.txt' is the viral dynamics model written in mlx format.
3. 'run.mlxtran' is the main monolix format that internally calls both the files above and runs NLME fitting.
4. Sub-folder 'run' is the actual monolix output which was used in the manuscript.
5. Once monolix converges wrt population distribution estimates of fit parameters, multiple individuals who are sampled from the distributions are simulated (in the sub-folder '\run\IndividualParameters'). 'Sel_Indiv_Params.csv' compiles the list of individuals whose parameters when input into the viral dynamics model best-captured individual clinical data from 'Neant_Kissler_Wolfel_input_monolix.csv' (Figure S7 in the paper).