% https://sysbiochalmers.github.io/Human-GEM-guide/gem_extraction/#prepare-the-transcriptomic-data
clear;
clc;
% load Human1 model
load('Human-GEM.mat')
load('C:\Users\prins\git\Human-GEM\model\prepData.mat')
% The second flag indicates if the model should be converted to gene symbols from ENSEMBL. This has to be decided at this point.
%prepData = prepHumanModelForftINIT(ihuman, false);
%save('prepData.mat', 'prepData')

%%

% read expression data and set 'data threshold'
data = read_expression_data("C:\Users\prins\git\\Human1_RPE_PR\RPE_PR data\test_data.xlsx")
data.threshold = 1;

for i = 1 : length(data.tissues)
    disp(string(data.tissues(i)))
    test_models_1step{i} = ftINIT(prepData, data.tissues{i}, [], [], data, {}, getHumanGEMINITSteps('1+0'), false, true);
end

save('models_1step_SysGO.mat', 'test_models_1step')

%%
baseModel = prepData.refModel;
models=models_1step

% now build a matrix saying which reactions are on
compMat = false(length(baseModel.rxns), length(models));

for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,models{i}.rxns);
end

% run t-sne
rng(1);  %set random seed to make reproducible
proj_coords = tsne(double(compMat.'), 'Distance', 'hamming', 'NumDimensions', 2, 'Exaggeration', 6, 'Perplexity', 3);

% export to R
d = struct();
d.tsneX = proj_coords(:, 1);
d.tsneY = proj_coords(:, 2);
save('TSNE.mat', 'd');


%%

gscatter(proj_coords(:,1),proj_coords(:,2),data.tissues')