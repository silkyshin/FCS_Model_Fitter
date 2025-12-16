A bunch of ODEs according to the oligomerization models I thought of fitting my FCS into.

Each folders in the ODE directory and their description:
forwardandreverse - fits both forward and reverse rates of each step
reduced - fits the forward rate and equilibrium constants at each step. Reverse rate is defined as k_fwd/k_eq

Will add other parameters/reduced versions in the future. You need to copy the models in the folders into the ODE directory and specify which one you want to run on ODE_solver.m

I might automate it later
