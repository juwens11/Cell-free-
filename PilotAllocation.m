function [allo_pilot_vec, produ_pilot_vec]= PilotAllocation(K,Vecteur_pilot,legth_tau)
   %%Random pilot allocation 
   allo_pilot_vec=randi([1 legth_tau], K, 1);
   produ_pilot_vec=[];
   for ii=1:K
       for jj=1:K
           if (allo_pilot_vec(ii)==allo_pilot_vec(jj))
           produ_pilot_vec(ii,jj)= 1 ;
           else
           produ_pilot_vec(ii,jj)=0 ;
           end 
       end     
   end

end