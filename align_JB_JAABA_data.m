function dat_JAABA=align_JB_JAABA_data(dat_JB,dat_JAABA)
%this function will eliminate some data from the JAABA array as it contains
%slightly more animal number than JB array
fn=fieldnames(dat_JAABA);
w=1;
for i=1:length(dat_JB.AN)
    while dat_JB.AN(i,1)~=dat_JAABA.AN{w,1} %find the aminal number that is in JAABA but not in JB and delete it 
        for j=1:length(fn)
            dat_JAABA.(fn{j})(w)=[];
        end
    end
    w=w+1;
end

if length(dat_JB.AN)<length(dat_JAABA.AN)
    idx=length(dat_JAABA.AN):-1:length(dat_JB.AN)+1;
    for j=1:length(fn)
        dat_JAABA.(fn{j})(idx)=[];
    end
end

end

