function [cost,extra] = netpc(x, Varmax, T, simtime,returnDetails)

    % Extract decision variables
    N_SMR   = x(1);
    N_solar = x(2);
    N_wind  = x(3);
    N_hydro = x(4);
%     penalty = x(5);
%     N_H2    = x(5);   
     BESSsize = x(5);
%     Extra = zeros(288,1);

    % Power generation per hour (288 points)
    P_solar = N_solar * T.PowerS;
    P_wind  = N_wind  * T.PowerW;
    P_hydro = N_hydro * T.PowerH;    % Already normalized to 1 MW base
    P_smr   = N_SMR   * T.PowerNuc;
    % The max powers of each unit
    
    P_solarmax = Varmax(2) * T.PowerS;
    P_windmax  = Varmax(3) * T.PowerW;
    P_hydromax = Varmax(4) * T.PowerH;    % Already normalized to 1 MW base
    P_smrmax   = Varmax(1) * T.PowerNuc;
    
    % At all times
    P_solarac = min(P_solar,P_solarmax);
    P_windac  = min(P_wind,P_windmax);
    P_hydroac = min(P_hydro,P_hydromax);
    P_smrac   = min(P_smr,P_smrmax);
    % Total supply per hour
    P_total = P_solarac + P_windac + P_hydroac + P_smrac;

     E_BESS = 0.9*BESSsize;
     BESSsocmin = 0.1*BESSsize;
     BESSsocmax = 0.9*BESSsize;
         
     P_totalf=P_total;

    LPSPF    = zeros(simtime,1);
    SEFF     = zeros(simtime,1);
    penaltyg = zeros(simtime,1);
     penaltyLP=zeros(simtime,1);
     penaltySF=zeros(simtime,1);
%      pensup = 0;
%     mismatch = zeros(simtime,1);
%     violation = zeros(simtime,1);

    for i = 2:length(T.P_el)
%      if P_totalf(i) < 1.2*T.P_el
%         penaltyd = penaltyd + (T.P_el(i) - P_totalf(i))^2 * 1e8;
%      end
    GRF_f = sum(P_totalf(i)) * 100/ sum(T.P_el(i)) ;

    end
    if (GRF_f < 105 || GRF_f > 125)  % if supply is more than 125% of demand
    penaltyg = penaltyg + (GRF_f - 100).^2 * 1e12;
    end 
    penalty = 0;
    c=0;
    d=0;
    BESSstate = zeros(288,1);
    for t = 2:length(T.P_el)
    demand = T.P_el(t);
    supply = P_totalf(t);
    
    if supply > demand
    % Surplus case ? Charge 
    surplus = supply - demand;
%     for s = 1:300
    charge = min(surplus * 0.85,BESSsocmax-E_BESS);
    E_BESS = E_BESS + charge;
%     c=charge;
%     end
    P_totalf(t) = P_totalf(t)-charge;
    end
    if supply < demand
    % Deficit case ? Discharge battery
    deficit = demand - supply;
    discharge = min(deficit/0.85,E_BESS-BESSsocmin);
    recovered = discharge * 0.85;
    E_BESS = E_BESS - discharge;
%     d = discharge;
    if E_BESS == BESSsocmin
       P_totalf(t)=P_totalf(t);
    else
    P_totalf(t) = supply + recovered; 
    end
%     end
%     Still deficit? Add penalty
    
    end
    BESSstate(t) = (E_BESS*100)/BESSsize;
    LPSPF = (sum(max(0, T.P_el(t) - P_totalf(t)+BESSsocmin))) * 100 / sum(T.P_el(t));
    SEFF  = (sum(max(0, P_totalf(t) - T.P_el(t)-BESSsocmax))) * 100 / sum(P_totalf(t));
%     if (P_totalf(t)<T.P_el(t) || P_totalf(t)>T.P_el(t))
%         penalty = penalty+1e8;
%     end
     if (P_totalf(t)> demand || P_totalf(t) < demand)
         penalty = penalty + 1e12;
     end
     Diff = P_totalf(t)-demand;
    end

    if (LPSPF > 5)
        penaltyLP = penaltyLP + (LPSPF - 5).^2 * 1e12;
    end
    if (SEFF > 10 || SEFF < 0)
        penaltySF = penaltySF+ (SEFF - 10).^2 * 1e12; % add penalty for violating LPSP and SEFF beyond limit
    end
%======Display of Ratios of Availability===================%
    disp(['Total GRF_f = ', num2str(GRF_f)]);
    disp(['Total LPSPF = ', num2str(LPSPF)]);
    disp(['Total SEFF = ', num2str(SEFF)]);
%==========================================================%
    % NPC calculation (simple example)
    % SMR, Solar, Wind, Hydro, Electrolyzer, FC, BESS
    cap_costs = [15, 1.2, 1.13, 2.5,0.398];  % $/unit (SMR, solar, wind, hydro/MW)
    CapMMR =0;
    for MMR=1:N_SMR
        CapMMR=CapMMR+(MMR)^(-0.2345)*cap_costs(1)*10000000;
    end
    NPC = CapMMR+N_solar*cap_costs(2)*600+N_wind*cap_costs(3)*500+N_hydro*cap_costs(4)*917235+...
    BESSsize*cap_costs(5);
   Repcost = N_solar*1*600*(1/(1+1.44)^25)...
             +N_wind*1.13*500*(1/(1+1.44)^25)...
             +BESSsize*0.398*(1/(1+1.44)^32);
   OandMcost = N_SMR*0.35*10000000*((1.06^40)-1)/(.06*(1.06)^40)...
              +N_solar*0.012*600*((1.06^40)-1)/(.06*(1.06)^40)...
              +N_wind*0.048*500*((1.06^40)-1)/(.06*(1.06)^40)...
              +N_hydro*0.1*917235*((1.06^40)-1)/(.06*(1.06)^40)...
              +BESSsize*0.01*((1.06^40)-1)/(.06*(1.06)^40);
   Fuelcost = (1*10^-5)*((1.06^40)-1)/(0.06*(1.06)^40)*10000000*300000;
   LTrep = 25*floor(40/25);
   LTrem = 25-(40-LTrep);
   Salvcost = (1*LTrem/25)+(1.13*LTrem/25)*10000000*N_SMR;
   Decomcost = (5*10^-6)*((1.06^40)-1)/(.06*(1.06^25))*10000000*300000;
   Refcost = N_SMR*(1/(1.06^44))*20000000*4;
   Difference = P_totalf-T.P_el;
   penaltybess = 0;
   if (BESSsize == 0)
       penaltybess = penaltybess + 1e8;
   end

    % Objective
    cost = NPC+Repcost+OandMcost+Fuelcost-Salvcost+Decomcost+Refcost+sum(penaltyg)...
        +sum(penaltyLP)+sum(penaltySF)+penaltybess+penalty;
    LCOE=cost-(sum(penaltyg)...
        +sum(penaltyLP)+sum(penaltySF)+penaltybess+penalty);
%   for LCOE = 1:length(T.P_el)
    COE = (LCOE/sum(P_totalf))*(0.07*(1.07)^40)/((1.07)^40-1)*1000;
%   end
    disp(['Total Cost= ', num2str(cost)]);
    disp(['SMR = ' num2str(N_SMR),' N_solar =' num2str(N_solar), ' N_wind = ' num2str(N_wind)...
    ,' N_hydro = ' num2str(N_hydro)]);
    disp([' BESSsize = ' num2str(BESSsize)]);
    disp([' E_BESS = ' num2str(E_BESS)]);
    disp(['Power output= ', num2str(sum(P_totalf)),'Load= ', num2str(sum(T.P_el))]);     
     disp(['penaltyd= ', num2str(sum(Difference))]);  
    disp(['COE=  ', num2str(COE)]);
    disp(['penaltyg= ', num2str(sum(penaltyg))]); 
    disp(['penaltyLP= ', num2str(sum(penaltyLP))]); 
    disp(['penalty = ', num2str(penalty)]);

    if nargin > 4 && returnDetails
        extra.P_totalf = P_totalf;
        extra.T.P_el = T.P_el;
        extra.GRF_f = GRF_f;
        extra.LPSPF = LPSPF;
        extra.SEFF = SEFF;
        extra.BESSstate = BESSstate;
        extra.penaltyg = penaltyg;
        extra.penalty = penalty;
        extra.Diff = Diff;
        % Add whatever else you need
    else
        extra = [];
    end
end


