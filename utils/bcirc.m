function BB = bcirc(AA)
    [mm,pp,nn] = size(AA);
    
    BB = [];
    
    for kk=1:nn
        BB = [BB; AA(:,:,kk)];
    end
    newCol = BB;
    for kk=2:nn
        newCol = circshift(newCol, mm);
        BB = [BB newCol];
    end
end