%%
clear all
load toycon1.mat; %load model
changeCobraSolver('gurobi'); %change solver (change to whatever solver you are using)
%%% TO MAKE FIGURES THE PACKAGE gramm and struct2csv must be added to the matlab path. The
%%% directory can be cloned from https://github.com/piermorel/gramm.git.
%%% Documetation is available with the command doc gramm. For struct2csv,
%%% documentation and download is available at
%%% https://www.mathworks.com/matlabcentral/fileexchange/34889-struct2csv.
%%% The foler must also be added to the path.
%%
toycon1_rxn_decompositon = struct(); %create structure
smat = full(model.S)  %grab smatrix
[numMet,numReact] = size(smat); %get number of reactions and metabolites
count =1 ; %counter for structure indexing
for j = 1:numReact
   k = find(smat(:,j)); %find indexes of non zero entries in a given column of smatrix (reaction)
   for l = 1:length(k)
        met_id(count) = model.mets(k(l));
        rxn_id(count) = model.rxns(j); %assign characteristics
        coeff(count) = smat(k(l),j);
        temp = char(met_id(count));
        met_compartment(count) = string(temp(5));  
        met_compound(count) = string(temp(1:3));
        cpd_name(count) = model.metNames(k(l));
        met_name(count) = strcat(model.metNames(k(l)), '[' , met_compartment(count) , ']'); 
        count = count + 1;
   end
end

%Create structure fields and observations
toycon1_rxn_decompositon.met_id = met_id;
toycon1_rxn_decompositon.rxn_id = rxn_id;
toycon1_rxn_decompositon.coeff = coeff;
toycon1_rxn_decompositon.met_compartment = met_compartment;
toycon1_rxn_decompositon.met_compount = met_compound;
toycon1_rxn_decompositon.cpd_name = cpd_name;
toycon1_rxn_decompositon.met_name = met_name;
disp(toycon1_rxn_decompositon)

%Output into text file
file = fopen('toycon1_rxn_decompositon.txt','w');
fprintf(file,'%s ',string(fieldnames(toycon1_rxn_decompositon)));
fprintf(file,'\n%s %s %s %s %s %s %s',[string(met_id);string(rxn_id);floor(coeff);string(met_compartment);string(met_compound)...
    ;string(cpd_name);string(met_name)]);
fclose(file);
%% Make S matrices txt files
%output S matrix with names
file = fopen('toycon1_smatrix_names.txt','w');
smatrix_names = [string(model.metNames),smat];
smatrix_names = [" ",string(model.rxnNames)';smatrix_names];
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',smatrix_names');
fclose(file);

%output S matrix with IDs
file = fopen('toycon1_smatrix_id.txt','w');
smatrix_id = [string(model.mets),smat];
smatrix_id = [" ",string(model.rxns)';smatrix_id];
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',smatrix_id');
%% Make rxn info struct
toycon1_rxn_info = struct();
%assign rxn_id,rxn_name,lower and upper bound fields
toycon1_rxn_info.rxn_id = model.rxns;
toycon1_rxn_info.rxn_name = model.rxnNames;
toycon1_rxn_info.lb = model.lb;
toycon1_rxn_info.ub = model.ub;
%construct reactions
for i = 1 :numReact
    sym = "-->";
    if model.rev(i) == 1%determine directionality symbol
        sym = "<==>";
    end
    temp = sym;
    k = find(smat(:,i)); 
    fcount = 0;
    rcount = 0;
    for j = 1 : length(k)
        comp = char(model.mets(k(j)));
        comp  = strcat('[',comp(5),']'); %get compartment
        if smat(k(j),i) < 0 %add onto reactants
            if fcount > 0
                temp = strcat(strcat(string(abs(smat(k(j),i))),model.metNames(k(j)),comp),strcat(' + ',temp));
            else
                temp = strcat(strcat(string(abs(smat(k(j),i))),model.metNames(k(j)),comp),temp);
            end
            fcount = fcount + 1;
        else  %add to products
            if rcount > 0
                temp = strcat(strcat(temp,' + '),strcat(string(abs(smat(k(j),i))),model.metNames(k(j))),comp);
            else
                temp = strcat(temp,strcat(string(abs(smat(k(j),i))),model.metNames(k(j))),comp);
            end
            rcount = rcount + 1;
        end
    end
    reactForm(i) = strrep(temp,'1',''); %remove 1 coefficent and add to vector
end
toycon1_rxn_info.rxn_formula = reactForm'; %assign reaction formula fields
disp(toycon1_rxn_info)

%% make S matrix pdf

pos = smat > 0;
temp = find(pos == 1)
ycoord = abs(numMet-mod(temp,numMet));
xcoord = floor(temp./numMet);
labels = smat(pos);
color = ones(length(xcoord),1);
neg = smat < 0;
temp = find(neg == 1);
xcoord = floor([xcoord ; temp./numMet]);
ycoord = mod([ycoord ; abs(numMet-mod(temp,numMet))],numMet);
labels = [labels ; smat(neg)];
color = [color;ones(length(temp),1)./2];


h = histogram2(xcoord,ycoord,'DisplayStyle','tile','ShowEmptyBins','on')