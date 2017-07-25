%%
clear all
load toycon1.mat; %load model
changeCobraSolver('gurobi'); %change solver (change to whatever solver you are using)
%%% TO MAKE FIGURES THE PACKAGE gramm must be added to the matlab path. The
%%% directory can be cloned from https://github.com/piermorel/gramm.git.
%%% Documetation is available with the command doc gramm. 
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
fig = figure
[Posrow,Poscol] = find(smat < 0 );%find negative coefficient matrix indices
ycoord = numMet - Posrow; %transform to xy coordinates
xcoord = numReact - Poscol;
labels = string(smat(smat < 0)); %make coefficent labels
labels = char(labels);
histogram2(xcoord,ycoord,'DisplayStyle','tile','ShowEmptyBins','off');%plot negative coefficient values in histogram
hold on;
text(xcoord,ycoord,labels,'color','white');%add labels
[Posrow,Poscol] = find(smat > 0);%repeat for positive coefficients
ycoord = numMet - Posrow;
xcoord = numReact - Poscol;
labels = char(string(smat(smat > 0)));
histogram2([xcoord;xcoord],[ycoord;ycoord],'DisplayStyle','tile','ShowEmptyBins','off');
text(xcoord,ycoord,labels,'color','white');
map = [1,1,0;0,0,1;1,0,0];%make color map
colormap(map);%apply color map
xticklabels(flip(model.rxnNames,1))%create xaxis tick labels
xtickangle(45); %rotate labels
for x = 1:numMet
    temp = char(model.mets{x}) ;
    metCompartments{x} = temp(length(temp)-2:length(temp));%create ylabels
end
yticklabels(strcat(flip(model.metNames),flip(metCompartments')))%add y tick labels
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

saveas(fig,'toycon1_smatrix','pdf')%save as pdf
%%      Perform flux variabilitiy analysis with percentage of objective function %%%%%%%%

fva_pct_result = ef_tbl_fva(0,model,toycon1_rxn_info,0);%perform initial fva
for i = 1:20
   fva_pct_result = ef_tbl_fva(i*5,model,fva_pct_result,1); %perform for all percentages
end

file = fopen('toycon1_fva_result_percentage.txt','w');%open file
fprintf(file,'rxn_id fva_lb fva_ub rxn_name lb ub rxn_formula fva_pct fva_req fva_on\n')%print headers
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',[string(fva_pct_result.rxn_id),string(fva_pct_result.fva_lb),...%print file
    string(fva_pct_result.fva_ub),string(fva_pct_result.rxn_name),string(fva_pct_result.lb),string(fva_pct_result.ub),...
    string(fva_pct_result.rxn_formula),string(fva_pct_result.fva_pct),string(fva_pct_result.fva_req),string(fva_pct_result.fva_on)]')

fclose(file);%close file

%%
%Plot percentage
fig = figure;
j = 1;
for rx = ["R1","R2"] %reactions to plot
    subplot(1,2,j)
    map = [.2,.2,.2;0,0,1;1,0,0];%make color map
    colormap(map);%apply color map
    k = strcmp(char(fva_pct_result.rxn_id) ,rx); %extract indexes for matching reactions
    xcoordl = fva_pct_result.fva_lb(k) ; %gather bounds for x coordinates
    xcoordu = fva_pct_result.fva_ub(k) ;
    ycoord = fva_pct_result.fva_pct(k); %gather pct for y coordinates
    color = coloring(fva_pct_result.fva_on(k),fva_pct_result.fva_req(k)); %make coloring based on fva_req and fva_off
    scatter(xcoordu,ycoord,[100],color,'filled','d') %Plot points
    hold on
    scatter(xcoordl,ycoord,[100],color,'filled','s')
    for i = 1:length(ycoord)
        switch color(i)
            case 1
                temp = 'red';
            case 0
                temp = 'black';
            case .5
                temp = 'cyan';
        end
        line([xcoordl(i),xcoordu(i)],[ycoord(i),ycoord(i)],'color',temp) %Draw connecting lines
    end
    j = j +1;
    xlabel('Units of Flux Through Reaction')
    ylabel('Required Flux Through Reaction')%label
    xlim([0,2])  %set xaxis limits
    temp = fva_pct_result.rxn_name(k);
    t = title(string(temp(1)));
    pos = get(t,'Position');
    pos(2) = pos(2) + 3;
    set(t,'Position',pos);
    grid on;
end
%save figure
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
saveas(fig,'toycon1_fva_percentage','pdf')

%%
% Peform FVA again but with ATP increments this time
fva_inc_result = ef_tbl_fva(0,model,toycon1_rxn_info,0);%perform initial fva
for i = 1:16
   fva_inc_result = ef_tbl_fva(100*i/16,model,fva_inc_result,1); %perform for all increments
end

file = fopen('toycon1_fva_result_increment.txt','w');%open file
fprintf(file,'rxn_id fva_lb fva_ub rxn_name lb ub rxn_formula fva_pct fva_req fva_on\n')%print headers
fprintf(file,'%s %s %s %s %s %s %s %s %s %s\n',[string(fva_inc_result.rxn_id),string(fva_inc_result.fva_lb),...%print file
    string(fva_inc_result.fva_ub),string(fva_inc_result.rxn_name),string(fva_inc_result.lb),string(fva_inc_result.ub),...
    string(fva_inc_result.rxn_formula),string(fva_inc_result.fva_pct),string(fva_inc_result.fva_req),string(fva_inc_result.fva_on)]')

fclose(file);%close file

%%
%Plot increment
fig = figure;
j = 1;
for rx = ["R1","R2"] %reactions to plot
    subplot(1,2,j)
    map = [.2,.2,.2;0,0,1;1,0,0];%make color map
    colormap(map);%apply color map
    k = strcmp(char(fva_inc_result.rxn_id) ,rx); %extract indexes for matching reactions
    xcoordl = fva_inc_result.fva_lb(k) ; %gather bounds for x coordinates
    xcoordu = fva_inc_result.fva_ub(k) ;
    ycoord = fva_inc_result.fva_pct(k).*32/100; %gather pct for y coordinates
    color = coloring(fva_inc_result.fva_on(k),fva_inc_result.fva_req(k)); %make coloring based on fva_req and fva_off
    scatter(xcoordu,ycoord,[100],color,'filled','d') %Plot points
    hold on
    scatter(xcoordl,ycoord,[100],color,'filled','s')
    for i = 1:length(ycoord)
        switch color(i)
            case 1
                temp = 'red';
            case 0
                temp = 'black';
            case .5
                temp = 'cyan';
        end
        line([xcoordl(i),xcoordu(i)],[ycoord(i),ycoord(i)],'color',temp) %Draw connecting lines
    end
    j = j +1;
    xlabel('Units of Flux Through Reaction')
    ylabel('Required # of ATP produced')%label
    xlim([0,2])  %set xaxis limits
    temp = fva_inc_result.rxn_name(k);
    t = title(string(temp(1)));
    pos = get(t,'Position');
    pos(2) = pos(2) + 1;
    set(t,'Position',pos);
    grid on;
end
%save figure
set(fig,'Units','Inches');
pos = get(fig,'Position'); %https://www.mathworks.com/matlabcentral/answers/12987-how-to-save-a-matlab-graphic-in-a-right-size-pdf
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]) 
saveas(fig,'toycon1_fva_increment','pdf')
%%%%%%%Perform Single Gene Deletions %%%%%%
%%
model = buildRxnGeneMat(model);
[~,~,~,~,~,sol] = singleGeneDeletion(model);
sol = sol(length(sol(:,1)),:)
model.genes'

