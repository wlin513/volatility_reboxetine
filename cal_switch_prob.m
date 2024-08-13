%calculate stay probability when Win and Loss appeared together with the left option or not at all respectively for each block
function out=cal_switch_prob(information,choice,resp_made)
    out=struct;

    sumofsame=0;
    sumofBoth=0;
    sumofNothing=0;
    sumofLossalone=0;
    sumofWinalone=0;
    sumofwin=0;
    sumofloss=0;
    sumofnowin=0;
    sumofnoloss=0;

    tempcount_same=0;
    tempcount_Both=0;
    tempcount_Nothing=0;
    tempcount_Winalone=0;
    tempcount_Lossalone=0;
    tempcount_win=0;
    tempcount_loss=0;
    tempcount_nowin=0;
    tempcount_noloss=0;
    tempcount_all=0;
    
    for i=1:length(information)-1
       if resp_made(i+1)
         if  choice(i)==1
            if information(i,1)==information(i,2)
                    sumofsame=sumofsame+1;
                if choice(i+1)==choice(i)
                    tempcount_same=tempcount_same+1;
                end
            end
            
            if information(i,1)/information(i,2)==1
                    sumofBoth=sumofBoth+1;
                if choice(i+1)==choice(i)
                    tempcount_Both=tempcount_Both+1;
                     tempcount_all=tempcount_all+1;
                end
            end
            
            if isnan(information(i,1)/information(i,2))
                    sumofNothing=sumofNothing+1;
                if choice(i+1)==choice(i)
                    tempcount_Nothing=tempcount_Nothing+1;
                      tempcount_all=tempcount_all+1;
                end
            end
            
            if isinf(information(i,1)/information(i,2))
                    sumofWinalone=sumofWinalone+1;
                if choice(i+1)==choice(i)
                    tempcount_Winalone=tempcount_Winalone+1;
                      tempcount_all=tempcount_all+1;
                end
            end
            
            if information(i,1)/information(i,2)==0
                    sumofLossalone=sumofLossalone+1;
                if choice(i+1)==choice(i)
                    tempcount_Lossalone=tempcount_Lossalone+1;
                      tempcount_all=tempcount_all+1;
                end
            end
            
            if information(i,1)==1
                    sumofwin=sumofwin+1;
                if choice(i+1)==choice(i)
                    tempcount_win=tempcount_win+1;
                end
            end
            
            if information(i,1)==0
                    sumofnowin=sumofnowin+1;
                if choice(i+1)==choice(i)
                    tempcount_nowin=tempcount_nowin+1;
                end
            end
         
            if information(i,2)==1
                    sumofloss=sumofloss+1;
                if choice(i+1)==choice(i)
                    tempcount_loss=tempcount_loss+1;
                end
            end  
            
            if information(i,2)==0
                    sumofnoloss=sumofnoloss+1;
                if choice(i+1)==choice(i)
                    tempcount_noloss=tempcount_noloss+1;
                end
            end   
         end
   
         
         
        if choice(i)==0
            if information(i,1)==information(i,2)
                    sumofsame=sumofsame+1;
                if choice(i+1)==choice(i)
                    tempcount_same=tempcount_same+1;
                end
            end
            
            if information(i,1)/information(i,2)==1
                    sumofNothing=sumofNothing+1;
                if choice(i+1)==choice(i)
                    tempcount_Nothing=tempcount_Nothing+1;
                    tempcount_all=tempcount_all+1;
                end
            end
            
            if isnan(information(i,1)/information(i,2))
                    sumofBoth=sumofBoth+1;
                if choice(i+1)==choice(i)
                    tempcount_Both=tempcount_Both+1;
                    tempcount_all=tempcount_all+1;
                end
            end
            
            if isinf(information(i,1)/information(i,2))
                    sumofLossalone=sumofLossalone+1;
                if choice(i+1)==choice(i)
                    tempcount_Lossalone=tempcount_Lossalone+1;
                    tempcount_all=tempcount_all+1;
                end
            end
            
            if information(i,1)/information(i,2)==0
                    sumofWinalone=sumofWinalone+1;
                if choice(i+1)==choice(i)
                    tempcount_Winalone=tempcount_Winalone+1;
                    tempcount_all=tempcount_all+1;
                end
            end
            
            if information(i,1)==0
                    sumofwin=sumofwin+1;
                if choice(i+1)==choice(i)
                    tempcount_win=tempcount_win+1;
                end
            end
            
            if information(i,1)==1
                    sumofnowin=sumofnowin+1;
                if choice(i+1)==choice(i)
                    tempcount_nowin=tempcount_nowin+1;
                end
            end
            
            if information(i,2)==0
                    sumofloss=sumofloss+1;
                if choice(i+1)==choice(i)
                    tempcount_loss=tempcount_loss+1;
                end
            end 
            
            if information(i,2)==1
                    sumofnoloss=sumofnoloss+1;
                if choice(i+1)==choice(i)
                    tempcount_noloss=tempcount_noloss+1;
                end
            end 
        end
      end
    end  
    
    staycount_same=tempcount_same;
    staycount_Both=tempcount_Both;
    staycount_Nothing=tempcount_Nothing;
    staycount_Winalone=tempcount_Winalone;
    staycount_Lossalone=tempcount_Lossalone;    
    staycount_win=tempcount_win;
    staycount_loss=tempcount_loss;  
    staycount_nowin=tempcount_nowin;
    staycount_noloss=tempcount_noloss;
    staycount_all=tempcount_all;

    out.switchprob_same=(1-staycount_same/sumofsame);
    out.switchprob_Both=(1-staycount_Both/sumofBoth);
    out.switchprob_Nothing=(1-staycount_Nothing/sumofNothing);
    out.switchprob_Winalone=(1-staycount_Winalone/sumofWinalone);
    out.switchprob_Lossalone=(1-staycount_Lossalone/sumofLossalone);    
    out.switchprob_win=(1-staycount_win/sumofwin);
    out.switchprob_loss=(1-staycount_loss/sumofloss); 
    out.switchprob_nowin=(1-staycount_nowin/sumofnowin);
    out.switchprob_noloss=(1-staycount_noloss/sumofnoloss);
    out.switchprob_all=(1-staycount_all/sum(resp_made));
    
   