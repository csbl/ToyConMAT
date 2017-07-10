clear all
load toycon1.mat;
changeCobraSolver('gurobi');
sol = optimizeCbModel(model);