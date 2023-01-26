function data = read_expression_data(filename,threshold)
    % read expression data from excel file with gene name in the first column and
    % expression data in the other columns (with header ids)

    if ~exist('threshold','var')
     % use 1 as default threshold
      threshold = 1;
    end
    t = readtable(filename, "UseExcel", false);
    % extract the tissue and gene names
    [~, n] = size(t);
    data.genes = t{:, 1}; % gene names
    data.tissues = t.Properties.VariableNames(3:n); % sample (tissue) names
    data.levels = t{:, 3:n}; % gene TPM values
end
