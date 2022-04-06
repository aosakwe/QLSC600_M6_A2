function dat_JAABA=arrange_JAABA_data(dat_JAABA)
%this function will simply concatenate the cells in each fields, so that it
%will have the same strcture as JB data
fn=fieldnames(dat_JAABA);
for i=1:length(fn)
    dat_JAABA.(fn{i})=vertcat(dat_JAABA.(fn{i}){:});% combine the cells 
for j=1:length(dat_JAABA.(fn{i}))
    dat_JAABA.(fn{i}){j,1}= dat_JAABA.(fn{i}){j,1}';
end

end
clear i j
end