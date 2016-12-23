%function [] = Reactor_plot(t_final, n_th_init, n_fa_init, V, m_tot, U5_pour, U8_pour, Pu9_pour, Poison_pour, T_PF, T_rt, Lambda_BC_thermal, Lambda_BC_fast)
function [] = GUI_1()
close all
    fig = figure('units','normalized','outerposition',[0.33 0.3 0.38 0.6],'name','Interface de la simulation de réacteur nucléaire 0D')
    text1 = uicontrol( fig , 'style' , 'text' , 'position' , [125,370,340,30] ,...
        'string' , 'Modélisation d''un réacteur 0D' , 'fontsize' , 15 )
    text2 = uicontrol( fig , 'style' , 'text' , 'position' , [50,320,500,30] ,...
        'string' , 'Veuillez donner une valeur aux paramètres suivants' , 'fontsize' , 15 )
    
    text3 = uicontrol ( fig , 'style' , ' edit' , 'position', [170,240,50,15] , 'Max' , 1 , 'string' , '100' );
    ui3 = uicontrol ( fig , 'style' , ' text' , 'position', [10,240,150,15] , 'string' , 'Durée de la simulation' );
    unit3 = uicontrol ( fig , 'style' , ' text' , 'position', [230,240,30,15] , 'string' , '[s]' );
    
    text4 = uicontrol ( fig , 'style' , ' edit' , 'position', [440,240,50,15] , 'Max' , 1 , 'string' , '25' );
    ui4 = uicontrol ( fig , 'style' , ' text' , 'position', [280,240,150,15] , 'string' , 'Masse totale de combustible' );
    unit4 = uicontrol ( fig , 'style' , ' text' , 'position', [500,240,30,15] , 'string' , '[kg]' );
    
    text5 = uicontrol ( fig , 'style' , ' edit' , 'position', [170,200,50,15] , 'Max' , 1 , 'string' , '1e15');
    ui5 = uicontrol ( fig , 'style' , ' text' , 'position', [10,175,150,45] , 'string' , 'Nombre initial de neutrons thermiques' );
    unit5 = uicontrol ( fig , 'style' , ' text' , 'position', [230,200,30,15] , 'string' , '[#]' );
    
    text6 = uicontrol ( fig , 'style' , ' edit' , 'position', [440,200,50,15] , 'Max' , 1 , 'string' , '0' )
    ui6 = uicontrol ( fig , 'style' , ' text' , 'position', [280,192,150,30] , 'string' , 'Nombre intitial de neutrons rapides' )
    unit6 = uicontrol ( fig , 'style' , ' text' , 'position', [500,200,30,15] , 'string' , '[#]' )
    
    text7 = uicontrol ( fig , 'style' , ' edit' , 'position', [170,147,50,20] , 'Max' , 1 , 'string' , '3' );
    ui7 = uicontrol ( fig , 'style' , ' text' , 'position', [10,142,150,30] , 'string' , 'Pourcentage massique d''U235 dans le combustible' );
    unit7 = uicontrol ( fig , 'style' , ' text' , 'position', [230,150,30,15] , 'string' , '[%]' );
    
    text8 = uicontrol ( fig , 'style' , ' edit' , 'position', [440,147,50,20] , 'Max' , 1 , 'string' , '97' );
    ui8 = uicontrol ( fig , 'style' , ' text' , 'position', [280,142 ,150,30] , 'string' , 'Pourcentage massique d''U238 dans le combustible' );
    unit8 = uicontrol ( fig , 'style' , ' text' , 'position', [500,150,30,15] , 'string' , '[%]' );
    
    text9 = uicontrol ( fig , 'style' , ' edit' , 'position', [170,108,50,15] , 'Max' , 1 , 'string' , '0' );
    ui9 = uicontrol ( fig , 'style' , ' text' , 'position', [10,100,150,30] , 'string' , 'Pourcentage massique d''U239 dans le combustible' );
    unit9 = uicontrol ( fig , 'style' , ' text' , 'position', [230,108,30,15] , 'string' , '[%]' );
    
    text10 = uicontrol ( fig , 'style' , ' edit' , 'position', [440,108,50,15] , 'Max' , 1 , 'string' , '5' );
    ui10 = uicontrol ( fig , 'style' , ' text' , 'position', [265,100,170,30] , 'string' , 'Pourcentage molaire de poison dans les produits de fission' )
    unit10 = uicontrol ( fig , 'style' , ' text' , 'position', [500,108,30,15] , 'string' , '[%]' );
        
    bp1 = uicontrol ( fig , 'style' , 'push' , 'position' , [420 10 60 30 ] ,...
        'string' , 'Démarrer' , 'callback' , @(bp1,eventdata) Reactor_model(str2double(get(text3,'String')),str2double(get(text5,'String')), str2double(get(text6,'String')), str2double(get(text4,'String')), str2double(get(text7,'String')), str2double(get(text8,'String')), str2double(get(text9,'String')), str2double(get(text10))))
    get(bp1,eventdata)
end
    

function []=GUI_2()

fig1=figure('units','normalized','outerposition',[0 0 1 1]);

text1 = uicontrol( fig1 , 'style' , 'text' , 'position' , [175,880,300,40] ,...
    'string' , 'Résulats' , 'fontsize' , 30 )

%Graphe T-s%
subplot ( 'Position' , [ .40 .6 .2 .3 ] ) ;
%Cloche%
T_0 = 1e-4; %[°C] car XSteam ne trouve rien pour T=0°C

S = linspace(0,XSteam('sV_T',T_0),200); %vecteur s dans la cloche
T = arrayfun(@(s) XSteam('Tsat_s',s), S);
plot(S,T,'k');
hold on;
plot(s_12,T_12,'blue')
plot(s_cBP,T_cBP,'blue')
plot(s_vBP,T_vBP,'blue')
plot(s_cHP,T_cHP,'blue')
plot(s_sBP,T_sBP,'blue')
plot(s_vHP,T_vHP,'blue')
plot(s_sHP,T_sHP,'blue')
plot(s_34,T_34,'blue')
plot(s_45,T_45,'blue')
plot(s_51,T_51,'blue')

text(etat1.s,etat1.t,'\leftarrowEtat 1~2')
text(etat3.s,etat3.t,'\leftarrowEtat 3')
text(etat4.s,etat4.t,'\leftarrowEtat 4')
title('Graphe T-s')
xlabel('Entropie [J/kgK]')
ylabel('Température [K]')

%Graphe h-s%
subplot ( 'Position' , [ .40 .13 .2 .3 ] ) ;
%dessine le diagramme (h,s) d'une centrale TV

%cloche TS
T_0 = 1e-2; %[°C]

S = linspace(XSteam('sL_T',T_0),XSteam('sV_T',T_0),200);
P = arrayfun(@(s) XSteam('psat_s',s), S);
H = arrayfun(@(p,s) XSteam('h_ps',p,s), P, S);

plot(S,H,'k');
hold on;
plot(s_12,h_12,'blue')
plot(s_cBP,h_cBP,'blue')
plot(s_vBP,h_vBP,'blue')
plot(s_cHP,h_cHP,'blue')
plot(s_sBP,h_sBP,'blue')
plot(s_vHP,h_vHP,'blue')
plot(s_sHP,h_sHP,'blue')
plot(s_34,h_34,'blue')
plot(s_45,h_45,'blue')
plot(s_51,h_51,'blue')

text(etat1.s,etat1.h,'\leftarrowEtat 1~2')
text(etat3.s,etat3.h,'\leftarrowEtat 3')
text(etat4.s,etat4.h,'\leftarrowEtat 4')
title('Graphe h-s')
xlabel('Entropie [J/kgK]')
ylabel('Enthalpie [kJ/kg]')


%%%%%%%%%%%
%Pie Chart%
%%%%%%%%%%%

%Energetique%

P=[double(P_eTG)*1000 double(P_eTV) double(P_fmec) double(P_cond_en) double(P_echap_en)]

label={sprintf('Puissance effective TG \n %0.1f MW ',P_eTG)...
    sprintf('Puissance effective TV \n %0.1f MW ',P_eTV/1e3)...
    sprintf('Pertes mécaniques \n %0.1f MW ',P_fmec/1e3)...
    sprintf('Pertes au condenseur: \n %0.1f MW ',P_cond_en/1e3)...
    sprintf('Pertes à l''échappement: \n %0.1f MW ',P_echap_en/1e3)};
subplot ( 'Position' , [ .66 .515 .3 .45 ] ) ;
pie(P,label);
title(sprintf('Distribution du flux énergétique avec puissance primaire de %0.1f  MW',P_prim_en/1e3 ));

%Exergetique%

P=[double(P_eTG)*1000 double(P_eTV)  double(P_fmec) double(P_cond_ex) double(P_rotex) double(P_echap_ex) double(P_irr_therm) double(P_irr_comb)]

label={sprintf('Puissance effective TG \n %0.1f MW ',P_eTG)...
    sprintf('Puissance effective TV \n %0.1f MW ',P_eTV/1e3)...
    sprintf('Pertes mécaniques \n %0.1f MW ',P_fmec/1e3)...
    sprintf('Pertes au condenseur: \n %0.1f MW ',P_cond_ex/1e3)...
    sprintf('Irréversibilités au complexe rotorique: \n %0.1f MW ',P_rotex/1e3)...
    sprintf('Pertes à l''échappement: \n %0.1f MW ',P_echap_ex/1e3)...
    sprintf('Irréversibilité du transfert thermique \n %0.1f MW ',P_irr_therm/1e3)...    
    sprintf('Irréversibilités de la combustion: \n %0.1f MW ',P_irr_comb/1e3)};

subplot ( 'Position' , [ .66 0.02 .3 .45 ] ) ;
pie(P,label);
title(sprintf('Distribution du flux exergétique avec puissance primaire de %0.1f  MW',P_prim_ex/1e3 ));

%%%%%%%%%%%%%%%%%%%%%%%
%Tableaux des résulats%
%%%%%%%%%%%%%%%%%%%%%9%
%Caractéristiques des états$
%Pour la TG
pg=[p1g;p2g;p3g;p4g;p5g];
tg=[t1g;t2g;t3g;t4g;t5g];
hg=[h1g;h2g;h3g;h4g;h5g];
sg=[s1g;s2g;s3g;s4g;s5g];
eg=[e1g;e2g;e3g;e4g;e5g];
Etatsg={'1g';'2g';'3g';'4g';'5g'};

Table = table(pg,tg,hg,sg,eg,'RowNames',Etatsg)


text1 = uicontrol( fig1 , 'style' , 'text' , 'position' , [170,686,330,40] ,...
    'string' , 'Caractéristiques des états de la TG' , 'fontsize' , 15 )
t1 = uitable(fig1);
t1.Data = table2cell(Table);
t1.Position = [110 561 425 132];
t1.ColumnName = {'p [kPa]','T [°C]','h [kJ/kg]','s [kJ/kgK]','e [kJ/kg]'};

%Pour la TV

p=[etat1.p;etat2.p;etat3.p;etat4.p;etat5.p];
t=[etat1.t;etat2.t;etat3.t;etat4.t;etat5.t];
x=[etat1.x;etat2.x;etat3.x;etat4.x;etat5.x];
h=[etat1.h;etat2.h;etat3.h;etat4.h;etat5.h];
s=[etat1.s;etat2.s;etat3.s;etat4.s;etat5.s];
e=[etat1.e;etat2.e;etat3.e;etat4.e;etat5.e];
Etats={'1';'2';'3';'4';'5'};

Table = table(p,t,x,h,s,e,'RowNames',Etats)


text1 = uicontrol( fig1 , 'style' , 'text' , 'position' , [170,500,330,40] ,...
    'string' , 'Caractéristiques des états de la TV' , 'fontsize' , 15 )
t1 = uitable(fig1);
t1.Data = table2cell(Table);
t1.Position = [110 375 425 132];
t1.ColumnName = {'p [bar]','T [°C]','x [-]','h [kJ/kg]','s [kJ/kgK]','e [kJ/kg]'};
% %Rendements$
% text2 = uicontrol( fig1 , 'style' , 'text' , 'position' , [170,300,300,40] ,...
%     'string' , 'Rendements' , 'fontsize' , 15 )
% ETAt=[double(eta_cyclen); double(eta_mec); double(eta_toten); double(eta_rotex); double(eta_cyclex); double(eta_combex); double(eta_totex)]
% ETA=table(ETAt,'RowNames',{'eta_cyclen'; 'eta_mec'; 'eta_toten'; 'eta_rotex'; 'eta_cyclex'; 'eta_combex'; 'eta_totex'})
% t2 = uitable(fig1);
% t2.Data = table2cell(ETA)%{eta_cyclen; eta_mec; eta_toten; eta_rotex; eta_cyclex; eta_combex; eta_totex}
% t2.Position = [210 355 200 160];
% t2.RowName = {'eta_cyclen','eta_mec','eta_toten','eta_rotex','eta_cyclex','eta_combex','eta_totex'};
% t2.ColumnName= {''}
% 
% %Puissances$
% text3 = uicontrol( fig1 , 'style' , 'text' , 'position' , [40,310,300,40] ,...
%     'string' , 'Puissances [MW]' , 'fontsize' , 15 )
% PUISSANCEt=[double(P_mT); double(P_mC); double(P_fmec); double(P_prim_en); double(P_echap_en); double(P_prim_ex); double(P_irr_comb);double(P_echap_ex);double(P_irr_tc)]/1000
% PUISSANCE=table(PUISSANCEt,'RowNames',{'P_mT'; 'P_mC'; 'P_fmec'; 'P_prim_en'; 'P_echap_en'; 'P_prim_ex'; 'P_irr_comb';'P_echap_ex';'P_irr_tc'})
% t3 = uitable(fig1);
% t3.Data = table2cell(PUISSANCE)%{eta_cyclen; eta_mec; eta_toten; eta_rotex; eta_cyclex; eta_combex; eta_totex}
% t3.Position = [80 115 205 195];
% t3.RowName = {'P_mT'; 'P_mC'; 'P_fmec'; 'P_prim_en'; 'P_echap_en'; 'P_prim_ex'; 'P_irr_comb';'P_echap_ex';'P_irr_tc'};
% t3.ColumnName= {''}

%Débits$
text4 = uicontrol( fig1 , 'style' , 'text' , 'position' , [310,320,300,40] ,...
    'string' , 'Débits [kg/s]' , 'fontsize' , 15 )
DEBITSt=[double(m_a); double(m_c); m3; m4]
DEBITS=table(DEBITSt,'RowNames',{'Air'; 'Carburant'; 'Vapeur au point 3'; 'Vapeur au point 4'})
t4 = uitable(fig1);
t4.Data = table2cell(DEBITS)%{eta_cyclen; eta_mec; eta_toten; eta_rotex; eta_cyclex; eta_combex; eta_totex}
t4.Position = [315 210 280 110];
t4.RowName = {'Air'; 'Carburant'; 'Vapeur au point 3'; 'Vapeur au point 4'};
t4.ColumnName= {''}

%Point de rosée%
text5 = uicontrol( fig1 , 'style' , 'text' , 'position' , [360,160,200,40] ,...
    'string' , 'Temp. de rosée [°C]' , 'fontsize' , 15 )
ROSEE=table(T_rosee,'RowNames',{'T_rosee'})
t5 = uitable(fig1);
t5.Data = table2cell(ROSEE);
t5.Position = [375 100 165 58];
t5.RowName = {'T_rosee'};
t5.ColumnName= {''}

%Boutons%
bp1 = uicontrol ( fig1 , 'style' , 'push' , 'position' , [70 50 200 30 ] ,...
    'string' , 'Nouvelle TGV' , 'callback' , @(bp1,eventdata)GUI_2(3))

bp2 = uicontrol ( fig1 , 'style' , 'push' , 'position' , [300 50 200 30 ] ,...
    'string' , 'Retour au menu principal' , 'callback' , @(bp2,eventdata)GUI())

end

