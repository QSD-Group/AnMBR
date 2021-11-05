tic

designs=fieldnames(output_E);
output_E.net_E=table();

for i=1:length(designs)
    output_E.(designs{i}).Energy = array2table(output_E.(designs{i}).Energy);
    output_E.net_E{i,1}=median(output_E.(designs{i}).Energy.Var3); %Calculates average and standard dev. and median from all trials
end
output_E.net_E.Properties.RowNames=designs;

toc