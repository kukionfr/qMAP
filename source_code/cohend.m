function d = cohend(matA,matB)
    meanA = mean(matA,1);meanB = mean(matB,1);
    stdA = std(matA,1);stdB = std(matB,1);
    d = abs(meanA-meanB)./((stdA+stdB)*0.5);
end