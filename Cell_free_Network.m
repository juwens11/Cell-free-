function [Pathloss_matrix ] =Cell_free_Network(Rc,M,K,d_1,d_0)


pentrationlossdB = 20; %Penetration loss (indoor users)
shadowFadingBS = 8; %Standard deviation (in dB) of shadow fading for all users to the macro BS
noiseFigure = 7; %dB
subcarrierBandwidth = 15e3;
%noiseFloordBm =  noiseFigure; %Noise floor in dBm

%% Pathloss parameter
Bandwidth = 20e6; % Bandwidth in Hz
%noiseFigure = 5;%dB
noiseFloordBm =-90; % -174+10*log10(Bandwidth) + noiseFigure;

vector_angle_users= (2*pi).*unifrnd(0,1,[K,1]);
vector_module_users= (Rc).*unifrnd(0,1,[K,1]);
Coordinate_matrix_users=zeros(K,2);


vector_angle_ap= (2*pi).*unifrnd(0,1,[M,1]);
vector_module_ap= (Rc).*unifrnd(0,1,[M,1]);
Coordinate_matrix_APs=zeros(M,2);

for i=1:M
    Coordinate_matrix_APs(i,1)=vector_module_ap(i)*cos(vector_angle_ap(i));
    Coordinate_matrix_APs(i,2)=vector_module_ap(i)*sin(vector_angle_ap(i));
end




for i=1:K
    Coordinate_matrix_users(i,1)=vector_module_users(i)*cos(vector_angle_users(i));
    Coordinate_matrix_users(i,2)=vector_module_users(i)*sin(vector_angle_users(i));
end





Distances_matrix=zeros(K,M) ;
for i=1:K
    for j=1:M
        Distances_matrix(i,j)=sqrt((Coordinate_matrix_users(i,1)-Coordinate_matrix_APs(j,1))^2+(Coordinate_matrix_users(i,2)-Coordinate_matrix_APs(j,2))^2);
    end
end
% figure
% 
% scatter(Coordinate_matrix_users(:,1),Coordinate_matrix_users(:,2))
% hold on
% scatter(Coordinate_matrix_APs(:,1),Coordinate_matrix_APs(:,2))
% leg1 = legend('Users','APs');
% set(leg1,'Location','Southeast')


Pathloss_matrix=zeros(K,M) ;

for i=1:K
    for j=1:M
        if   Distances_matrix(i,j)>d_1
            Pathloss_matrix_dB=  -35.7-35*log10(Distances_matrix(i,j))+ shadowFadingBS *randn(1,1)- noiseFloordBm;
            
        elseif Distances_matrix(i,j)< d_0
            
            Pathloss_matrix_dB= -81.2; %  15*log10(d_1)+20*log10(d_0)  -59.5+shadowFadingSCA_outcluster*randn(1,1)+ noiseFloordBm;
        else
            
            Pathloss_matrix_dB=-61.2-20*log10(Distances_matrix(i,j));
            
        end
        Pathloss_matrix(i,j)=10^((Pathloss_matrix_dB)/10);
        
    end
end

end 