function [Debit, PILOT]= power_alloca(ASSign_pilo, Debit_group,Nouveau_W,C_total,Mean_op, Rho_p,Tau_p,Const,Rho_d,Trans_ante,Large_Scale_Fading_matrix,K,M) 


Debit_Achievable_group=Debit_group;

Cost=7;


 



% PILOT=zeros(K,1);
% ite_p=0;
% X=Mean_op*ones(C_total,1);
% for group=1:C_total    
%     AA=find(X==group);
%     vec_el=[];
%     for jj=1:length(AA)
%         if(ite_p < Tau_p)
%             ite_p=ite_p+1;
%             PILOT(AA(jj))=ite_p;
%         else
%             PILOT(AA(jj))=0;
%         end 
%     end 
% end 
% YY=find(PILOT~=0);
% for k=1:length(YY)
%     Interf=zeros(Tau_p,1);
%     for pil=1:Tau_p
%           for jj=1:K
%                         if (jj~= YY(k)) %(ind_user~= Y_user(k))
%                             if (ASSign_pilo(jj)==pil)
%                                 Interf(pil)= Interf(pil) +   Nouveau_W(jj,YY(k)) ; %sum(Large_Scale_Fading_matrix(k,:)); % +  Nouveau_W(ind_user,Y_user(k));
%                             else
%                                 Interf(pil)= Interf(pil);
%                             end 
%                             
%                         end
%         end
%     end
%        [valmi, indmi]=min(Interf);
%        PILOT(YY(k))= indmi; % PILOT(INDEX_I(indmi));    
% end 
%    
%   
%              
% 
%         NOUvew_numkgroup=zeros(K,M, C_total);
%         for group=1:C_total
%             for m=1:M
%                 for k=1:K
%                     if (Mean_op(k,group)==1)
%                         Deno_sum=0;
%                         for jj=1:K
%                             if (jj~= k)
%                                 if(PILOT(jj)==PILOT(k))
%                                     Deno_sum = Deno_sum +  Large_Scale_Fading_matrix(jj,m);
%                                 end
%                             end
%                         end
%                        NOUvew_numkgroup(k,m,group) = Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)^2/(Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)+Tau_p*Rho_p*Deno_sum +1);
%                     end
%                 end
%             end
%         end
% 
%         Interference_jkgroup=zeros(K,K, C_total);
%         Debit_Achievable_group=zeros(K,C_total);
%         for group=1:C_total
%             for k=1:K
%                 if (Mean_op(k,group)==1)
%                     First_term=0;
%                     second_term=0;
%                     for jj=1:K
%                         if (jj~= k)
%                             if(Mean_op(jj,group)==1)
%                                 if(PILOT(jj)==PILOT(k))
%                                    First_term=Large_Scale_Fading_matrix(k,:)*diag(1./Large_Scale_Fading_matrix(jj,:))*(NOUvew_numkgroup(jj,:,group)');
%                                 end 
%                               % First_term=Large_Scale_Fading_matrix(k,:)*diag(1./Large_Scale_Fading_matrix(jj,:))*(nu_mkgroup(jj,:,group)');
%                                 second_term=Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(jj,:,group)');
%                                 
%                                 Interference_jkgroup(k,jj,group)=First_term + second_term;
%                             end
%                         end
%                     end
%                     Debit_Achievable_group(k,group)= Const*log2(1+ ( Rho_d*(Trans_ante*sum(NOUvew_numkgroup(k,:,group)))^2  )/(Rho_d*Trans_ante*Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(k,:,group)')+...
%                         1+ Rho_d*Trans_ante*sum(Interference_jkgroup(k,:,group))));      
%                 end  
%             end
%         end 
%   ASSign_pilo= PILOT;    
 
    TTT=0;
while(TTT < 60)
    TTT=TTT+1;
  Debit_Achievable_group(Debit_Achievable_group==0)=100;
   [val_row, indsrow]=min(Debit_Achievable_group);
   [Valcol, inde_col]=min(val_row);
    K_star=indsrow(inde_col);

 pil_assig=ASSign_pilo;
 Rece_inter_pilot_k=zeros(Tau_p,1);
  for pil=1:Tau_p
      pil_assig(K_star)=pil;
      nu_mk_pilotalloc=zeros(K,M);
      for k=1:K
          if (k~=K_star)
            for m=1:M
                        Deno_sum=0;
                        for jj=1:K
                            if (jj~=k)
                                if(pil_assig(jj)==pil_assig(k))
                                    Deno_sum = Deno_sum +  Large_Scale_Fading_matrix(jj,m);
                                end 
                            end
                        end
                        nu_mk_pilotalloc(k,m) = Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)^2/( Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m) +Tau_p*Rho_p*Deno_sum +1);
            end
          end
      end
       for jj=1:K
                        if (jj~= K_star)
                            if(pil_assig(jj)==pil)
                                Rece_inter_pilot_k(pil)=Rece_inter_pilot_k(pil)+ nu_mk_pilotalloc(jj,:)*Large_Scale_Fading_matrix(k,:)';
                               else
                                Rece_inter_pilot_k(pil)= Rece_inter_pilot_k(pil);
                            end 
                            
                        end 
       end 
  
  end
   
   
   
[Val_pro, ind_pilot]=min(Rece_inter_pilot_k);
ASSign_pilo(K_star)=ind_pilot;
 
 PILOT=ASSign_pilo;
 NOUvew_numkgroup=zeros(K,M, C_total);
        for group=1:C_total
            for m=1:M
                for k=1:K
                    if (Mean_op(k,group)==1)
                        Deno_sum=0;
                        for jj=1:K
                            if (jj~= k)
                                if(PILOT(jj)==PILOT(k))
                                    Deno_sum = Deno_sum +  Large_Scale_Fading_matrix(jj,m);
                                end
                            end
                        end
                       NOUvew_numkgroup(k,m,group) = Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)^2/(Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)+Tau_p*Rho_p*Deno_sum +1);
                    end
                end
            end
        end

        Interference_jkgroup=zeros(K,K, C_total);
        Debit_Achievable_group=zeros(K,C_total);
        for group=1:C_total
            for k=1:K
                if (Mean_op(k,group)==1)
                    First_term=0;
                    second_term=0;
                    for jj=1:K
                        if (jj~= k)
                            if(Mean_op(jj,group)==1)
                                if(PILOT(jj)==PILOT(k))
                                   First_term=Large_Scale_Fading_matrix(k,:)*diag(1./Large_Scale_Fading_matrix(jj,:))*(NOUvew_numkgroup(jj,:,group)');
                                end 
                              % First_term=Large_Scale_Fading_matrix(k,:)*diag(1./Large_Scale_Fading_matrix(jj,:))*(nu_mkgroup(jj,:,group)');
                                second_term=Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(jj,:,group)');
                                
                                Interference_jkgroup(k,jj,group)=First_term + second_term;
                            end
                        end
                    end
                    Debit_Achievable_group(k,group)= Const*log2(1+ (Cost* Rho_d*(Trans_ante*sum(NOUvew_numkgroup(k,:,group)))^2  )/(Rho_d*Trans_ante*Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(k,:,group)')+...
                        1+ Rho_d*Trans_ante*sum(Interference_jkgroup(k,:,group))));      
                end  
            end
        end 






end     
        
        
        
        
        
        
        
   Debit=Debit_Achievable_group;
end 