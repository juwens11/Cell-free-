
clear all; 
close all; 
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rc=500;
M=100;
K=30;
%% Pathloss parameter
d_1= 50;%m
d_0= 10; %m
Shstd=8;

ValeurRho_d=[0.2 0.4 0.6 0.8 1 1.2 1.4 1.6];
% for k=1:100
%    aphla(k)= 2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rho_d=1;
Rho_p=0.2;
Tau_c=200;
Tau_p=20;
Nombre_frequence=10;
Trans_ante=10;
BarRate=0.1;
Bandwidth=10;
Const1=1.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_user_group= 6;
C_total=5; 
Alpha_gt=5;
legth_tau=Tau_p;

[Vecteur_pilot]=pilotvec(legth_tau);

Num_AP=3:7; 
Const=1-Tau_p/Tau_c;
%%pilot vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOTAL=100;
Corel_SANSgROUP=[];
Corelation_AvecG=[];
Debit_with_group=[];
Debit_Sans_grouping=zeros(K,TOTAL);
Debit_Allocation=zeros(K,TOTAL);
Allocation_Bande=zeros(C_total, TOTAL);
MonteCARLOTAverage_down_proposed=zeros(length(Num_AP),1);
MONTERAverag_pro_randomFrequency=zeros(length(Num_AP),1);
Path_los=[];



Scheduling_Propo=[];
Lagragian_exing=[];
Power_SCA_freq=[];

mont=0;
while (mont <TOTAL)
    mont=mont+1;
    flagg=false;
       MM11=100;
       [Alloca_pilot, Norm_Pilot_product]=PilotAllocation(K,Vecteur_pilot,legth_tau);

       [Large_Scale_Fading_matrix ]= Cell_free_Network(Rc,MM11,K,d_1,d_0); %dB
       M=100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Corelation matrix without group
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W_int=zeros(K,K);
        Nouveau_W=zeros(K,K);
        for k=1:K
            Term1=0;
            Term2=0;
                    for m=1:M
                        Inter=0;
                        Inter=(Large_Scale_Fading_matrix(k,m)*Large_Scale_Fading_matrix(:,m))./(Large_Scale_Fading_matrix(k,m)+Large_Scale_Fading_matrix(:,m)+1);
                        Term1=Term1+ (Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)*(Large_Scale_Fading_matrix(:,m)))./(Tau_p*Large_Scale_Fading_matrix(k,m)+Tau_p*Rho_p*Large_Scale_Fading_matrix(:,m)+1);
                        Term2= Term2+ (Tau_p*Rho_p*Large_Scale_Fading_matrix(k,m)*(Large_Scale_Fading_matrix(:,m).^2))./(Tau_p*Large_Scale_Fading_matrix(k,m)+Tau_p*Rho_p*Large_Scale_Fading_matrix(:,m)+1);
                        W_int(k,:)=W_int(k,:) + Inter';
                    end
                    Term1=Trans_ante*Term1.^2;
                    Nouveau_W(k,:)=Term1'+ Term2';
        end
        Max_des2W=zeros(K,K);
       for k=1:K
            for jj=1:K
                if (jj~=k)
                        Max_des2W(k,jj)=max(Nouveau_W(k,jj), Nouveau_W(jj,k));
                end
            end
       end
        W_int=W_int-diag(diag(W_int));
        W_matrice=W_int;
        %W_matrice=W_matrice/(10^10);
        Tilde_W=[W_matrice zeros(K,1); zeros(1,K+1)];
        
        Nouveau_W=0.5*(Nouveau_W+Nouveau_W');
        Nouveau_W= Nouveau_W-diag(diag(Nouveau_W));
        
         Gros_Augment_W= zeros(C_total*K);
         Aug_max_W= zeros(C_total*K);
         NOU_aug_w= zeros(C_total*K);

        for k=1:K
            for j=1:C_total
                Gros_Augment_W((1+(j-1)*K): (K+ (j-1)*K),(1+(j-1)*K): (K+ (j-1)*K))= W_matrice;%+10^(5)*eye(K,K) ;
                Aug_max_W((1+(j-1)*K): (K+ (j-1)*K),(1+(j-1)*K): (K+ (j-1)*K))= Max_des2W;
                NOU_aug_w((1+(j-1)*K): (K+ (j-1)*K),(1+(j-1)*K): (K+ (j-1)*K))=Nouveau_W;
            end
        end
        
        
         Delta_k=zeros(K,C_total*K);
        for jj=1:K
            Theta=zeros(1,K);
            Theta(1,jj)=1;
            Delta_k(jj,:)=repmat(Theta,1,C_total);
        end
        
       vect_int=ones(1,C_total*K);
       Feasible_sol=[]; 
        
       
       Nouvelle_cout=[];
       Maxiii_object=[];
       Fonction_cout=[];
       for itera=1:1000; 
        Condit_exit=0;
        while (Condit_exit==0)
            Condit_exit=1;
             X_inter=randi([0 1], C_total*K, 1);
                     for k=1:K
                         if(Delta_k(k,:)*X_inter  > Alpha_gt)
                          Condit_exit=0;
                         end 
                     if (Delta_k(k,:)*X_inter == 0)
                       Condit_exit=0;
                     end 
                     end
                        for j=1:C_total
                            if ( sum(X_inter((1+(j-1)*K): (K+ (j-1)*K))) > Tau_p)
                           Condit_exit=0;
                            end 
                   end
        end 
           va_f=(X_inter)'*Gros_Augment_W*X_inter; 
        Fonction_cout=[Fonction_cout va_f];
        Nouvelle_cout=[Nouvelle_cout (X_inter)'*NOU_aug_w*X_inter];
        Feasible_sol=[Feasible_sol X_inter];
       end 
 
         
         [Len, Kon]=min(Nouvelle_cout);
         Mean_op=reshape(Feasible_sol(:,Kon),[K,C_total]);
         
        %% Cardinal of each group
        Mew_Card=sum(Mean_op,1);

        
Sol_Lagrangian_perit=zeros(length(ValeurRho_d),1);
Propo_sol_perit=zeros(length(ValeurRho_d),1);
Group_sca_frequency=zeros(length(ValeurRho_d),1);
for it_rho=1:length(ValeurRho_d)
 Rho_d=ValeurRho_d(it_rho);    
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%     Lagrangian method 
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Random feasible pilot allocation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Feasi_pilot=rand(K,Nombre_frequence);
Feasi_solution=Feasi_pilot;
Const=1-Tau_p/Tau_c;
Varsigma=50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update of the lagrangian multipliers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


Gamma_k=10*rand(K, 1);
Upsilonk=5*rand(K,1);
Lagr_A=zeros(K,Nombre_frequence);
Lagrangian_B=zeros(K,Nombre_frequence);
Lagrang_C=zeros(K,1);
for k=1:K
    for freq=1:Nombre_frequence
        Lagr_A(k,freq)=Gamma_k(k)*(BarRate-Throuphput(k,freq));
        Lagrangian_B(k,freq)=Feasi_pilot(k,freq)*Throuphput(k,freq)- Varsigma*(Feasi_pilot(k,freq)-Feasi_pilot(k,freq)^2);
    end 
    Lagrang_C=Upsilonk(k)*(1-sum(Feasi_pilot(k,:)));
end 
Lagrangian_funct=sum(sum(Lagr_A))-sum(sum(Lagrangian_B))+ sum(Lagrang_C);

Lagran_val_int=0;
val=1.13;
while (abs((Lagrangian_funct-Lagran_val_int)/Lagran_val_int) > 10^(-3))
    Lagran_val_int=Lagrangian_funct;
    inter_gammak=Gamma_k;
    interVALupsil=Upsilonk;

  
    
    [Derirate_j_k, NU_MKL,Throuphput]=DERIVATIVE_Rjl_xkl(Alloca_pilot,Rho_d, Const, Feasi_pilot,Rho_p,Tau_p,K,M,Large_Scale_Fading_matrix(:,1:M),Trans_ante,Nombre_frequence);
 
    
    Solution_lagran=10*ones(K,Nombre_frequence);
    for freq=1:Nombre_frequence
        for k=1:K
            AAA=0;
            for jj=1:K
                if (jj~= k)
                    AAA=AAA-(inter_gammak(jj)+1)*Derirate_j_k(jj,k,freq);
                else
                    AAA=AAA-(inter_gammak(k)+1)*Throuphput(k,freq);
                end 
            end 
           Solution_lagran(k,freq)=Feasi_pilot(k,freq)-0.05*(-interVALupsil(k)+ AAA +Varsigma*(1-2*Feasi_pilot(k,freq))) ;
           Solution_lagran(k,freq)=max(0,Solution_lagran(k,freq));
           Solution_lagran(k,freq)=min(1,Solution_lagran(k,freq));
        end 
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%    Compute the rate associated to this solution
    
    SINR_kl=zeros(K, Nombre_frequence);
    Achievable_rate=zeros(K, Nombre_frequence);

    for freq=1:Nombre_frequence
            for k=1:K
                First_term=0;
                for jj=1:K
                    if (jj~= k)
                       First_term= First_term + IKjl_k(k,jj)*Solution_lagran(jj,freq);
                    end 
                end 
                   SINR_kl(k,freq)= 1+ ( Ak(k) )/(  Bk(k)+  Rho_d*Trans_ante*First_term );
                   Achievable_rate(k,freq)=Const*log2(1+SINR_kl(k,freq)*val)/(10^2);
            end 
    end 
    %%  End of  achievable rate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute of tthe lagrangian 
    Lagr_A=zeros(K,Nombre_frequence);
    Lagrangian_B=zeros(K,Nombre_frequence);
    Lagrang_C=zeros(K,1);
    for k=1:K
        for freq=1:Nombre_frequence
            Lagr_A(k,freq)=Gamma_k(k)*(BarRate-Achievable_rate(k,freq));
            Lagrangian_B(k,freq)=Solution_lagran(k,freq)*Achievable_rate(k,freq)- Varsigma*(Solution_lagran(k,freq)-Solution_lagran(k,freq)^2);
        end
        Lagrang_C=Upsilonk(k)*(1-sum(Solution_lagran(k,:)));
    end
    Lagrangian_funct=sum(sum(Lagr_A))-sum(sum(Lagrangian_B))+ sum(Lagrang_C);
    %% END of the lagrangian 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Upsilonk=zeros(K,1);
    Gamma_k=zeros(K,1);
    for k=1:K
        %for l=1:Nombre_frequence
        Gamma_k(k)=max(0, inter_gammak(k)+0.5*(BarRate- Solution_lagran(k,:)*(Achievable_rate(k,:)')));
        % end
        Upsilonk(k)=max(0, interVALupsil(k)+ 0.5*(1-sum(Solution_lagran(k,:))));
    end
    Feasi_pilot=Solution_lagran;
end 


 
 
 
 
 
 
        
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Rate for each user in each group
        
      % [norm_new_alloca] = New_pilot_allocation(Group_Card, C_total, K);
       NOUvew_numkgroup=zeros(K,M, C_total);
        for group=1:C_total
            for m=1:M
                for k=1:K
                    if (Mean_op(k,group)==1)
                        Deno_sum=0;
                        for jj=1:K
                            if (jj~= k)
                                if(Alloca_pilot(jj)==Alloca_pilot(k))
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
       Rate_Achievable_group=zeros(K,C_total);
        for group=1:C_total
            for k=1:K
                if (Mean_op(k,group)==1)
                    First_term=0;
                    second_term=0;
                    for jj=1:K
                        if (jj~= k)
                            if(Mean_op(jj,group)==1)
                                if(Alloca_pilot(jj)==Alloca_pilot(k))
                                   First_term=Large_Scale_Fading_matrix(k,1:M)*diag(1./Large_Scale_Fading_matrix(jj,1:M))*(NOUvew_numkgroup(jj,1:M,group)');
                                end 
                                second_term=Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(jj,:,group)');
                                
                                Interference_jkgroup(k,jj,group)=First_term + second_term;
                            end
                        end
                    end
                    Rate_Achievable_group(k,group)= Const*log2(1+ ( Rho_d*(Trans_ante*sum(NOUvew_numkgroup(k,:,group)))^2  )/(Rho_d*Trans_ante*Large_Scale_Fading_matrix(k,1:M)*(NOUvew_numkgroup(k,:,group)')+...
                        1+ Rho_d*Trans_ante*sum(Interference_jkgroup(k,:,group))));      
                end  
            end
        end 
  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%  Successive convex approximation


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Successive convex approximation 
THETA=rand(Nombre_frequence,C_total, 10);
Optimal_obj=[];
for  con=1:10
Theta_intermediaire=rand(Nombre_frequence,C_total);

Rate_group=sum(Rate_Achievable_group,1);

Valeur_objective=sum(sum(Theta_intermediaire*Rate_group'))-100*sum( sum(Theta_intermediaire-Theta_intermediaire.^2      ));
Val_interm=0;
while (abs((Valeur_objective-Val_interm)/Val_interm)> 10^(-4))
    Val_interm=Valeur_objective;
cvx_begin quiet 
    variable  thetagroup(Nombre_frequence,C_total) 
    maximize(  sum(thetagroup*Rate_group')-100*sum( sum(thetagroup-Theta_intermediaire.^2 -2*Theta_intermediaire.*(thetagroup-Theta_intermediaire))))
    subject to 
          0 <= thetagroup  ;
           thetagroup <= 1 ;
          1 <= sum(thetagroup,1); 
        %  sum(thetagroup,2)<=1;
          for k=1:K 
               BarRate <= sum(thetagroup*(Rate_Achievable_group(k,:)'))  ;
          end 
cvx_end 
    thetagroup(thetagroup < 0)=0;
 Theta_intermediaire=thetagroup;
Valeur_objective=cvx_optval ;
if(cvx_optval==-Inf)
    flagg=1; 
    break
end 
end 
if(flagg)
    break
end 
Optimal_obj(con)=cvx_optval ;
THETA(:,:,con)=thetagroup;
end 
if(flagg)
    break
end 

[val_ob, ind_obj]=max(Optimal_obj);
OptimTheta=THETA(:,:,ind_obj);    
        
        
        
         Final_thetagroup=round(OptimTheta);
        Prop_debit_peruser=zeros(K,1);
        for k=1:K
           Prop_debit_peruser(k)=Const1*(1/Nombre_frequence)*sum((Final_thetagroup*(Rate_Achievable_group(k,:)'/(10^2))));
        end 
[Throughput, pi_assign]=power_alloca(Alloca_pilot, Rate_Achievable_group,Nouveau_W,C_total,Mean_op, Rho_p,Tau_p,Const,Rho_d,Trans_ante,Large_Scale_Fading_matrix,K,M) ;
Final_rate=zeros(K,1);
       for k=1:K
           Final_rate(k)=(1/Nombre_frequence)*sum((Final_thetagroup*(Throughput(k,:)')))/(10^2);
        end 


Sol_Lagrangian_perit(it_rho)=(sum(mean(Achievable_rate,2))+sum(Prop_debit_peruser))/2;
Propo_sol_perit(it_rho)=sum(Prop_debit_peruser); 
Group_sca_frequency(it_rho)=sum(Final_rate); 
  
   
         
          
         
          
         
         
     


 end 
     if(flagg)
            mont=mont-1;
          continue
     end
       

Scheduling_Propo=[ Scheduling_Propo Propo_sol_perit];
Lagragian_exing=[Lagragian_exing Sol_Lagrangian_perit];
Power_SCA_freq=[Power_SCA_freq Group_sca_frequency];
% save Scheduling_Lagragian_intermediare
end



grid on
hold on
plot(ValeurRho_d, 20*mean(Power_SCA_freq,2),'-v',ValeurRho_d, 20*mean(Scheduling_Propo,2),'-o',ValeurRho_d,20*mean(Lagragian_exing,2),'-d')
xlabel('Downlink Transmit Power (Watt)')
ylabel('Sum Downlink rate (Mbits/s)')
hold on
legend('Proposed SDR + SCA +Power ','Proposed SDR + SCA ','Lagrangian-based solution')
grid on

%save Lagrangian_scheduling_de_pilot_plus_scheduling