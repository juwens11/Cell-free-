function [Deri_rate_j_k, NU_mkl,Throuphput]=DERIVATIVE_Rjl_xkl(Alloca_pilot,Rho_d, Const, Feasi_solution,Rho_p,Tau_p,K,M,Large_Scale_Fading_matrix,Trans_ante,Nombre_frequence)

Constant=Const/log(2);
NU_mkl=zeros(K,M);
    for m=1:M
        for k=1:K
            Deno_sum=0;
            for jj=1:K
                if (jj~= k)
                    if(Alloca_pilot(jj)==Alloca_pilot(k))
                       Deno_sum = Deno_sum +  Large_Scale_Fading_matrix(jj,m);
                    end 
                end 
            end 
            NU_mkl(k,m) = Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)^2/(Rho_p*Large_Scale_Fading_matrix(k,m) +Tau_p*Rho_p*Deno_sum +1);
        end 
    end 
    Ak=zeros(K,1);
    Bk=zeros(K,1);
     for k=1:K
         Ak(k)=Rho_d*(Trans_ante*sum(NU_mkl(k,:)))^2 ;
         Bk(k)=Rho_d*Trans_ante*Large_Scale_Fading_matrix(k,1:M)*(NU_mkl(k,:)')  +1 ;
     end 

     
    IKjl_k=zeros(K,K);
        for k=1:K
            for jj=1:K
                if (jj~=k)
                    if(Alloca_pilot(jj)==Alloca_pilot(k))
                    Premier_term=Trans_ante*(Large_Scale_Fading_matrix(k,1:M)*diag(1./Large_Scale_Fading_matrix(jj,1:M))*(NU_mkl(jj,:)'))^2;
                    IKjl_k(k,jj)=Premier_term+ Large_Scale_Fading_matrix(k,1:M)*(NU_mkl(jj,:)');
                    else 
                     IKjl_k(k,jj)= Large_Scale_Fading_matrix(k,1:M)*(NU_mkl(jj,:)');
                    end 
                end
            end 
        end  
        
        

SINR_kl=zeros(K, Nombre_frequence);
Throuphput=zeros(K, Nombre_frequence);
Deno_SINRkl=zeros(K, Nombre_frequence);
for freq=1:Nombre_frequence
        for k=1:K
            First_term=0;
            for jj=1:K
                if (jj~= k)
                   First_term= First_term + IKjl_k(k,jj)*Feasi_solution(jj,freq);
                end 
            end 
               SINR_kl(k,freq)= 1+ ( Ak(k) )/(  Bk(k)+  Rho_d*Trans_ante*First_term );
               Deno_SINRkl(k,freq)= Bk(k)+  Rho_d*Trans_ante*First_term  ;
               Throuphput(k,freq)=Const*log2(1+SINR_kl(k,freq));
        end 
end 
Deri_rate_j_k=zeros(K,K, Nombre_frequence);
for freq=1:Nombre_frequence
        for k=1:K
            for jj=1:K
                if (jj~= k)
                   Deri_rate_j_k(k,jj,freq)=((Constant)*(-1)*(Ak(k)* IKjl_k(k,jj))/(Deno_SINRkl(k,freq)^2))/SINR_kl(k,freq);
                end 
            end 
        end 
end 


end 