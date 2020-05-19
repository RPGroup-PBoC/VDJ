for j=1:8
    lac=lacANAL_RMS{j};
    for i=1:size(lac,1)
        close all
        plot(lac(i,:))
        aa=[];
        aa=ginput(1);
            i
            aa
        lac(i,round(aa):end)=0;
    end
    lacANAL_RMS{j}=lac;
end