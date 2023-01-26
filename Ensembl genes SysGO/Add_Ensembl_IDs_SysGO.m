
% SysGO database 
SysGO = readtable('SysGO_Ensembl_expression_HUMAN1');

% reference table
genes = readtable('genes.xlsx');

% match gene names based on geneSymbols
for i = 1 : length(SysGO.GeneSymbol)
  idx = find(ismember(genes.geneSymbols,  SysGO.GeneSymbol(i)));
      if ~isempty(idx)
      SysGO.Ensembl(i) = genes.genes(idx);
      else % if no match found look for UNIprotID
          if  ~isempty(char(SysGO.uniprot_ids(i)))
           idx = find(ismember(genes.geneUniProtID,  SysGO.uniprot_ids(i)));
               if ~isempty(idx);
                SysGO.Ensembl(i) = genes.genes(idx);
               else
               end
          else
        end
    end
end

% could not match following Ensemble IDs
no_match = setdiff(genes.genes,SysGO.Ensembl(~cellfun(@isempty,SysGO.Ensembl)));

N_matches = sum(~cellfun(@isempty,SysGO.Ensembl));
disp( ['Could match ',num2str(N_matches),...
    ' Ensemble genes out of ', num2str(length(genes.genes))]);

N_double_Uniprot_matches = N_matches+length(no_match)-length(genes.genes);
disp( ['There were ',num2str(N_double_Uniprot_matches),...
    ' double matches based on UniProtID']);

% SysGO table with matched Ensembl IDs
writetable(SysGO,'SysGO_Ensembl.xlsx'); 

% table with unmatched Ensembl IDs
idx=[];
for i = 1 : length(no_match)
    idx(i) = find(ismember(genes.genes,no_match(i)));
end

unmatched_Ensembl = genes(idx,:)
writetable(unmatched_Ensembl,'unmatched_Ensembl.xlsx')
%%


% list all the alternative gene names and see whether there is overlap in
% the last few unmatched genes

% for unmatched genes, make list of gene aliases and gene symbols in Ensembl table
g=[];
for i = 1 : length(unmatched_Ensembl.genes)
     if  ~isempty(char(unmatched_Ensembl.geneAliases(i)))
     gi = split(unmatched_Ensembl.geneAliases(i),';');
     g  = vertcat(g,gi);
     end
     if  ~isempty(char(unmatched_Ensembl.geneSymbols(i)))
     gi = split(unmatched_Ensembl.geneSymbols(i),';');
     g  = vertcat(g,gi);
     end
end

% for all genes, make list of gene aliases and gene symbols in Ensembl table
g_sysgo=[];
for i = 1 : size(SysGO,1)
     if  ~isempty(char(SysGO.alias_symbol(i)))
     gi = split(SysGO.alias_symbol(i),'|');
     g_sysgo  = vertcat(g_sysgo,gi);
     end
     if  ~isempty(char(SysGO.prev_symbol(i)))
     gi = split(SysGO.prev_symbol(i),'|');
     g_sysgo  = vertcat(g_sysgo,gi);
     end
     if  ~isempty(char(SysGO.GeneSymbol(i)))
     gi = SysGO.GeneSymbol(i);
     g_sysgo  = vertcat(g_sysgo,gi);
     end
end

idx=[]
count=1;
for i = 1 : length(g)
      idx_i = find(ismember(g_sysgo,  g(i)));
      if ~isempty(idx);
      idx(count)=idx_i;
      end
      count+1;
end

