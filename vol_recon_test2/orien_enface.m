function ori2D = orien_enface(O3D,binsize)


         ori2D = zeros(size(O3D,1),size(O3D,2));
            for i=1:size(O3D,1)
                for j = 1:size(O3D,2)
                    yO=squeeze(O3D(i,j,:));
                    yO(yO==0)=[];
                    yO(isnan(yO))=[];

                    if yO
                        vectorO=-90:binsize:90;
                        nO=hist(yO,vectorO);
    %                     figure,hist(yO,vectorO);
                        [maxO,indexO]=max(nO);
                                               
                        ori2D(i,j)=vectorO(indexO);
                    end
                    clear yO
                    
                end
            end
            
end