%  Data Preparation for Random Forest and GAM
%  Johns Hopkins University
%  Shengzhi Luan
%  02.08.2023
% =========================================================================
Threshold = 0.9;
%  ------------------------------------------------------------------------
Feature_Table = readtable('Structure-Property-Data.csv', ...
                          'PreserveVariableNames',true);
Feature_Matrix = readmatrix('Structure-Property-Data.csv');
fprintf('-- Calculating the correlation between features ...\n')
Correlation = corr(Feature_Matrix(:,1:end-1));
[Pair_I,Pair_J] = find(abs(Correlation)>=Threshold & abs(Correlation)~=1);
Pair = unique(sort([Pair_I,Pair_J],2),'rows');
%  ------------------------------------------------------------------------
for i = 1:1:size(Pair,1)
    Text_Message = ['---- Correlated features found: ', ...
                    Feature_Table.Properties.VariableNames{Pair(i,1)}, ...
                    ' & ', ...
                    Feature_Table.Properties.VariableNames{Pair(i,2)}, ...
                    ' | Pearson correlation value as ', ...
                    num2str(Correlation(Pair(i,1),Pair(i,2)))];
    fprintf(strcat(Text_Message,'...\n'));
end
fprintf('-- First feature of each correlated pair is removed ...\n')
%  ------------------------------------------------------------------------
fprintf('-- Generating new feature matrix files ...\n')
mkdir('Refined Data');
cd('Refined Data');
Feature_Matrix(:,Pair(:,1)) = [];
Feature_Table(:,Pair(:,1)) = [];
writetable(Feature_Table,'Structure-Property-Refined-Data.csv');
Feature_Table_1P = Feature_Table([1:5:end],:);
writetable(Feature_Table_1P,'Structure-Property-Refined-1P-Data.csv');
cd ..
clear all
% =========================================================================