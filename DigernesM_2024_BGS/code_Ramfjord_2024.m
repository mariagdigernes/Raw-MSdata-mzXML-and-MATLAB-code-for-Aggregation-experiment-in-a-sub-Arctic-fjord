%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% csv assign version 7
%
% Originally written by Jeffrey Hawkes and modified by Maria G. Digernes        
%
% Feel free to distribute, but please acknowledge the programmers and curators!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Upload Ramfjord data
% define and load files
clear all

%reassign=0; % choose 0 to use old data, 1 to apply changes in raw data, total number of samples or sample order 

myFolder='/Users/mariagdigernes/OneDrive - NTNU/MATLAB/DigernesM_2024_BGS/RAWmzXMLfiles';
HostFolder='/Users/mariagdigernes/OneDrive - NTNU/MATLAB/DigernesM_2024_BGS';

cd(HostFolder) %command to go to folder with all files 
data_type=2;% 1=csv, 2=mzXML, 3=xlsx, 4=txt
ESI_mode=-1; %negative mode
[smpdat,names_samps] = xlsread('samplenamesRAMFJORD.xlsx');% creates two new variables smpdat (values only) and namessamps (names only).
nosamps=size(names_samps,1); %number of samples. calculates the number of samples on samplename file and stores on column 1
noreps=1;
sampnames=names_samps(:,1); % creates samp names as variable for long sample name (located on column 1 of samplenames file)
sampcodes=names_samps(:,2);  % creates samp codes as variable for short sample name (located on column 2 of samplenames file)

%sampcodes={'ESFA';'PLFA';'SRFA'}  
cal_peak=[];
cal_peak=369.11911; %remove if you do not need a linear internal cal
linear_cal_ppm=10;
fine_cal_ppm=3;
assign_ppm=0.7;
calchoice=0;

min_mass=150;
max_mass=800;

%elements: [min,max,exact mass,valence] this creates a table with each
%value specified below for each column
%(numer,:) means make the matrix on the right 
elements(1,:)=[4,50,12,4]; %C
elements(2,:)=[4,100,1.007825,1]; %H
elements(3,:)=[2,40,15.9949146,2]; %O
elements(4,:)=[0,2,14.003074,3]; %N
elements(5,:)=[0,1,31.972071,4]; %S
elements(6,:)=[0,0,30.973762,5]; %P
elements(7,:)=[0,0,79.916521,6]; %Se
elements(8,:)=[0,0,22.989770,1]; %Na
elements(9,:)=[0,0,34.96885269,7]; %Cl
elements(10,:)=[0,1,13.00335,4]; %C13

no_els=size(elements,1); %counts total number of elements and puts the number on column 1

%if reassign==1
assign %runs script called assign 
%else
%load('matlab') %runs all of the variables saved on matlab file (.mat)
%end

[smpdat,names_samps] = xlsread('samplenamesRAMFJORD.xlsx'); %added these next 5 lines so that it loads samplename file again when
%changed
nosamps=size(names_samps,1);%number of samples. calculates the number of samples on samplename file and stores on column 1
noreps=1;
sampnames=names_samps(:,1);
sampcodes=names_samps(:,2);  
samplenamesfilter=9; %column to plot or current working column from samplenames file

fontsizeall=16;
tannin=formula(:,3)>0.65&formula(:,2)<1.5; %JAH

%% Mass spectra
FigMS = figure(); %
FigMS.Color=([1 1 1]);
FigMS.Units='centimeters';
FigMS.Position=([5 5 100 100]); 
%
sample_startval_ms = 1 ;%added start sample number 
for i=1:16  %position to graph
        rep=1;
        subplot(8,2,i) %(rows,columns,position) number of plots 
        smp= i + sample_startval_ms - 1 ;%changed added to connect sample start val to smp
        datplot=avg_intensities(:,smp);
        stem(formula(:,1),datplot,'Marker','none') %formula contains all intensities for all 10 elements 
        %axis([150 850 0 prctile(datplot(datplot>0),99.9)])
        axis([150 850 0 7000])
        ylabel('ion abundance'),xlabel('m/z')
        title(sampcodes{smp});      
end 
%% Figure 7 mass spectra overlay
fontsizeall=14;

%September
FigMSOverlay = figure('Name','September Mass spectra');%
FigMSOverlay.Color=([1 1 1]);
FigMSOverlay.Units='centimeters';
FigMSOverlay.Position=([10 5 30 10]);

smp_under = 3; %smpdat row number
smp_over = 4; %last one you plot comes on top

datplot_under = avg_intensities(:,smp_under); %time 0
datplot_over = avg_intensities(:,smp_over); %time 1
stem(formula(:,1),datplot_under,'Marker','none','Color',"#FFC20A",'LineWidth',2); %orange t0
hold on
stem(formula(:,1),datplot_over,'Marker','none','Color',"#556B2F",'LineWidth',2); %green time 1
%axis([150 850 0 prctile(datplot_under(datplot_under>0),99.9)]) %xmin xmax ymin ymax
axis([150 850 0 6000 ]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')
legend('t0','t1');
title('September');
annotation('textbox',[.05 .5 .5 .5],'String','b)','EdgeColor','none', 'FontSize', 19) %(x y width height)
set(gca,'Fontsize',fontsizeall);

%December
FigMSOverlay = figure('Name','December Mass spectra');%
FigMSOverlay.Color=([1 1 1]);
FigMSOverlay.Units='centimeters';
FigMSOverlay.Position=([10 5 30 10]);

smp_under = 19; %smpdat row number
smp_over = 20; %last one you plot comes on top

datplot_under = avg_intensities(:,smp_under); %time 0
datplot_over = avg_intensities(:,smp_over); %time 1
stem(formula(:,1),datplot_under,'Marker','none','Color',"#FFC20A",'LineWidth',2); %orange t0
hold on
stem(formula(:,1),datplot_over,'Marker','none','Color',"#556B2F",'LineWidth',2); %green time 1
%axis([150 850 0 prctile(datplot_under(datplot_under>0),99.9)]) %xmin xmax ymin ymax
axis([150 850 0 6000 ]) %xmin xmax ymin ymax
ylabel('ion abundance'),xlabel('m/z')
legend('t0','t1');
title('December');
annotation('textbox',[.05 .5 .5 .5],'String','c)','EdgeColor','none', 'FontSize', 19) %(x y width height)
set(gca,'Fontsize',fontsizeall);




%% Figure S8 VK overlay 
%September
gridtransparency=0.3;

if true  
    
    FigVK = figure('Name','VK:September'); 
    FigVK.Color=([1 1 1]);
    FigVK.Units='centimeters';
    FigVK.Position=([5 3 20 15]);
    cut=0;
    ax=gca();
    grid(ax, 'on') ;
    set (ax, 'Fontsize',15);
    ax.LineWidth = 1.5;
    ax.GridColor = 'k';
    ax.GridAlpha = gridtransparency;
    sample_startval_vk = 1 ;
    
    for i= 1:4 %position on subplot. also total number of samples
        smp = i + sample_startval_vk -1; %added to connect smp to sample startval
        if (smpdat(i,1) ~= 0)% samples equal to 0 on smpdat are removed.  col 7 are 3 cruises on 1 plot.
            datplot=avg_intensities(:,smp); %
            sizex=datplot/max(datplot)*800; % 
            %
            a=datplot>cut; %all avg intensities from selected samples
            b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
            c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % N>0 S=0, CHON
            d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS
            
            formulatype= a; 
            hold on;
            if smpdat(smp,1) == 1  ||  smpdat(smp,1) ==2 ||  smpdat(smp,1) ==3 %time 0
                time0=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#FFC20A", "MarkerEdgeColor","#FFC20A", 'DisplayName','t0'); %yellow Sep
            elseif smpdat(smp,1) == 4 %time 1
                time1=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#556B2F","MarkerEdgeColor","#556B2F","MarkerFaceAlpha",0.5, "MarkerEdgeAlpha",0.5,'DisplayName','t1');%green sep
            end
            
            annotation('textbox',[.05 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 19) %(x y width height)
            axis([0 1 0 2.5])
            ylabel('H/C')
            %xlabel('O/C')
            title('September'); 
            %legend([time0,time1],'t0','t1','Location', 'northeast')
            %
        end
    end
    clear b c d cut formulatype

%October
   FigVK = figure('Name','VK:October'); 
    FigVK.Color=([1 1 1]);
    FigVK.Units='centimeters';
    FigVK.Position=([5 3 20 15]);
    cut=0;
    ax=gca();
    grid(ax, 'on') ;
    set (ax, 'Fontsize',15);
    ax.LineWidth = 1.5;
    ax.GridColor = 'k';
    ax.GridAlpha = gridtransparency;
    sample_startval_vk = 9 ;
    
    for i= 1:4 %position on subplot. also total number of samples
        smp = i + sample_startval_vk -1; %added to connect smp to sample startval
        if (smpdat(i,1) ~= 0)% t.
            datplot=avg_intensities(:,smp); %
            sizex=datplot/max(datplot)*800; % 
            %
            a=datplot>cut; %all avg intensities from selected samples
            b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
            c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % N>0 S=0, CHON
            d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS
            
            formulatype= a; 
            hold on;
            if smpdat(smp,1) == 1  ||  smpdat(smp,1) ==2 ||  smpdat(smp,1) ==3 %time 0
                time0=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#FFC20A", "MarkerEdgeColor","#FFC20A", 'DisplayName','t0'); %yellow Sep
            elseif smpdat(smp,1) == 4 %time 1
                time1=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#001f3F","MarkerEdgeColor","#001f3F","MarkerFaceAlpha",0.5, "MarkerEdgeAlpha",0.5,'DisplayName','t1'); %blue winter
            end
            
            annotation('textbox',[.05 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 19) %(x y width height)
            axis([0 1 0 2.5])
            %ylabel('H/C')
            %xlabel('O/C')
            title('October'); 
            %
            legend([time0,time1],'t0','t1','Location', 'northeast')
            %
        end
    end
    clear b c d cut formulatype


%December 
    FigVK = figure('Name','VK:December'); 
    FigVK.Color=([1 1 1]);
    FigVK.Units='centimeters';
    FigVK.Position=([5 3 20 15]);
    cut=0;
    ax=gca();
    grid(ax, 'on') ;
    ax.LineWidth = 1.5;
    ax.GridColor = 'k';
    ax.GridAlpha = gridtransparency;
    set (ax, 'Fontsize',15);
    sample_startval_vk = 17 ;
    
    for i= 1:4 
        smp = i + sample_startval_vk -1; 
        if (smpdat(i,1) ~= 0)% 
            datplot=avg_intensities(:,smp); 
            sizex=datplot/max(datplot)*800; 
            %
            a=datplot>cut; %all avg intensities from selected samples
            b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
            c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % N>0 S=0, CHON
            d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS
            
            formulatype= a;
            hold on;
            if smpdat(smp,1) == 1  ||  smpdat(smp,1) ==2 ||  smpdat(smp,1) ==3 %time 0
                time0=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#FFC20A", "MarkerEdgeColor","#FFC20A", 'DisplayName','t0'); %yellow Sep
            elseif smpdat(smp,1) == 4 %time 1
                time1=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#001f3F","MarkerEdgeColor","#001f3F","MarkerFaceAlpha",0.5, "MarkerEdgeAlpha",0.5,'DisplayName','t1'); %blue winter
            end
            
            annotation('textbox',[.05 .5 .5 .5],'String','(c)','EdgeColor','none', 'FontSize', 19) %(x y width height)
            axis([0 1 0 2.5])
            ylabel('H/C')
            xlabel('O/C')
            title('December'); 
            %legend([time0,time1],'t0','t1','Location', 'northeast')
            %
        end
    end
    clear b c d cut formulatype



%February 
    FigVK = figure('Name','VK:February'); 
    FigVK.Color=([1 1 1]);
    FigVK.Units='centimeters';
    FigVK.Position=([5 3 20 15]);
    cut=0;
    ax=gca();
    grid(ax, 'on') ;
    ax.LineWidth = 1.5;
    ax.GridColor = 'k';
    ax.GridAlpha = gridtransparency;
    set (ax,'Fontsize',15)
    sample_startval_vk = 25 ;
    
    for i= 1:4 %position on subplot. also total number of samples
        smp = i + sample_startval_vk -1; %
        if (smpdat(i,1) ~= 0)% samples equal to 0 on smpdat are removed. 
            datplot=avg_intensities(:,smp); %
            sizex=datplot/max(datplot)*800; 
            %
            a=datplot>cut; %all avg intensities from selected samples
            b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
            c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % N>0 S=0, CHON
            d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS
            
            formulatype= a;
            hold on;   
            if smpdat(smp,1) == 1  ||  smpdat(smp,1) ==2 ||  smpdat(smp,1) ==3 %time 0
                time0=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#FFC20A", "MarkerEdgeColor","#FFC20A", 'DisplayName','t0'); %yellow Sep
            elseif smpdat(smp,1) == 4 %time 1
                time1=scatter(formula(formulatype,3),formula(formulatype,2),sizex(formulatype),"MarkerFaceColor","#001f3F","MarkerEdgeColor","#001f3F","MarkerFaceAlpha",0.5, "MarkerEdgeAlpha",0.5,'DisplayName','t1'); %blue winter
            end
            
            annotation('textbox',[.05 .5 .5 .5],'String','(d)','EdgeColor','none', 'FontSize', 19) %(x y width height)
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('February'); 
            %legend([time0,time1],'t0','t1','Location', 'northeast')
            
        end
    end
    clear b c d cut formulatype

end 
%legend works if remove legend then plot then add legend after.



%% Figure S6 vK molecular weight

LOWMW=formula(:,1)<250;  
HIGHMW=formula(:,1)>570; 

%September low MW
FigLMW = figure('Name','September m/z filter');
FigLMW.Color=([1 1 1]);
FigLMW.Units='centimeters';
FigLMW.Position=([5 3 20 15]);
iter=1;
cut=0;
sample_first = 1 ;
sample_last = 3 ;
    for i= sample_first:sample_last
        plotposition = smpdat(i,1); %go to column 2 on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample ie. 1.2.3...!
            datplot=avg_intensities(:,i);
            sizex=datplot/max(datplot)*1000; %intensities of peaks affects dot size
            a=datplot>cut;
            b=datplot>cut&LOWMW==1; 
            hold on
            scatter(formula(a,3),formula(a,2),sizex(a),"MarkerFaceColor","#D3D3D3","MarkerEdgeColor","#D3D3D3") 
            scatter(formula(b,3),formula(b,2),sizex(b),"MarkerFaceColor","#007500","MarkerEdgeColor","#007500")
            axis([0 1 0 2.5])
            ylabel('H/C')
            xlabel('O/C')
            title('September');
            legend('t0','m/z < 250')
            set(gca,'FontSize', 19);
            annotation('textbox',[.015 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end
    
%December high MW 
FigHMW = figure('Name','December m/z filter');
FigHMW.Color=([1 1 1]);
FigHMW.Units='centimeters';
FigHMW.Position=([5 3 20 15]);
iter=1;
cut=0;
sample_first = 17 ;
sample_last = 19 ;
    for i= sample_first:sample_last
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample ie. 1.2.3...!
            datplot=avg_intensities(:,i);
            sizex=datplot/max(datplot)*1000; 
            a=datplot>cut;
            c=datplot>cut&HIGHMW==1; %HMW
            hold on
            scatter(formula(a,3),formula(a,2),sizex(a),"MarkerFaceColor","#D3D3D3","MarkerEdgeColor","#D3D3D3") 
            scatter(formula(c,3),formula(c,2),sizex(c),"MarkerFaceColor","#007500","MarkerEdgeColor","#007500")
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('December');
            legend('t0','m/z > 570')
            set(gca,'FontSize', 19);
            annotation('textbox',[.015 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end
clear b c d cut

%% Figure S7 vK t-peaks

% 184 t-peak formulas from Medeiros et al 2016

%Take chemical formula "CH4" and store as C=1 and H=4

%Formula in the format: CxHyNzOa

%Column 1: C
%Column 2: H
%Column 3: N
%Column 4: O
%Column 5: O/C
%Column 6: H/C
%Column 7: Molecular mass

[NUM,TXT,RAW]=xlsread('Tpeaks.xlsx');
Tpeaks_molecularformula = string(RAW(2:185,1));
Tpeaks_molecularmass = NUM(:,1);
Tpeaks_AI = NUM(:,2);
Tpeaks_OC = NUM(:,3);
Tpeaks_HC = NUM(:,4);

Tpeaks_table = zeros(numel(Tpeaks_molecularformula),7);
Tpeaks_table(:,5) = Tpeaks_OC;
Tpeaks_table(:,6) = Tpeaks_HC;
Tpeaks_table(:,7) = Tpeaks_molecularmass;


for i = 1:numel(Tpeaks_molecularformula)
    linedata = Tpeaks_molecularformula(i);
    linemass = Tpeaks_molecularmass(i);

    if contains(linedata,"O")
       tmp = split(linedata,"O"); % tmp(1) is remainder, tmp(2) is # of O
       linedata = tmp(1);
       O = str2num(tmp(2));
    else
       O = 0;
    end
    
    if contains(linedata,"N")
       tmp = split(linedata,"N"); % tmp(1) is remainder, tmp(2) is # of N
       linedata = tmp(1);
       N = str2num(tmp(2));
    else
        N = 0;
    end
    
    if contains(linedata,"H")
       tmp = split(linedata,"H"); % tmp(1) is remainder, tmp(2) is # of H
       linedata = tmp(1);
       H = str2num(tmp(2));
    else
        H = 0;
    end
    
    if contains(linedata,"C")
       tmp = split(linedata,"C"); % tmp(1) is remainder, tmp(2) is # of C
       linedata = tmp(1);
       C = str2num(tmp(2));
    else
        C = 0;
    end
       
    Tpeaks_table(i,1) = C;
    Tpeaks_table(i,2) = H;
    Tpeaks_table(i,3) = N;
    Tpeaks_table(i,4) = O;

end
 
elements(1,:)=[4,50,12,4]; %C
elements(2,:)=[4,100,1.007825,1]; %H
elements(3,:)=[2,40,15.9949146,2]; %O
elements(4,:)=[0,2,14.003074,3]; %N
elements(5,:)=[0,1,31.972071,4]; %S
elements(6,:)=[0,0,30.973762,5]; %P
elements(7,:)=[0,0,79.916521,6]; %Se
elements(8,:)=[0,0,22.989770,1]; %Na
elements(9,:)=[0,0,34.96885269,7]; %Cl
elements(10,:)=[0,1,13.00335,4]; %C13


% Make formulas into string from 2229 peaks (in all samples) and store in molecular file_complete 
molecularformula_complete = strings(size(elements_used,1),1);
for i = 1:size(molecularformula_complete,1)
    %elements_used=[C(:) H(:) O(:) N(:) S(:) P(:)];
   if (elements_used(i,1)+elements_used(i,10)) > 0 %C and C13 combined
       molecularformula_complete(i) = strcat(molecularformula_complete(i),"C",num2str((elements_used(i,1)+elements_used(i,10))));
   end  
   if elements_used(i,2) > 0
       molecularformula_complete(i) = strcat(molecularformula_complete(i),"H",num2str(elements_used(i,2)));
   end  
   if elements_used(i,3) > 0
       molecularformula_complete(i) = strcat(molecularformula_complete(i),"O",num2str(elements_used(i,3)));
   end  
   if elements_used(i,4) > 0
       molecularformula_complete(i) = strcat(molecularformula_complete(i),"N",num2str(elements_used(i,4)));
   end  
   if elements_used(i,5) > 0
       molecularformula_complete(i) = strcat(molecularformula_complete(i),"S",num2str(elements_used(i,5)));
   end  
end

%check if formula match from t-peaks in all samples
Tpeaks_formula_filter = false(size(molecularformula_complete,1),1);
for i = 1:size(Tpeaks_formula_filter,1)
   if max(strcmp(molecularformula_complete(i),Tpeaks_molecularformula))
       Tpeaks_formula_filter(i) = 1;
   end
end

Tpeaks_match = false(size(Tpeaks_molecularformula,1),1);
for i = 1:size(Tpeaks_match,1)
    if max(strcmp(molecularformula_complete,Tpeaks_molecularformula(i)))
        Tpeaks_match(i) = 1;
    end
end
    
%make vK figures 

% t-peaks vK diagram Sept t0
FigTpeaksSEP = figure('Name','September t-peaks');
subplot(4,2,1);
FigTpeaksSEP.Color=([1 1 1]);
FigTpeaksSEP.Units='centimeters';
FigTpeaksSEP.Position=([5 3 50 25]);
dotsize = 50;
cut= 0;
sample_first_SEP_t0 = 1 ; 
sample_last_SEP_t0 = 3;

    for i= sample_first_SEP_t0:sample_last_SEP_t0
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            ylabel('H/C'),xlabel('O/C')
            title('September t0');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t0',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.07 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end


% t-peaks vK diagram sep t1
subplot(4,2,2)
FigTpeaksFEB.Color=([1 1 1]);
dotsize = 50;
cut= 0;
sample_first_SEP_t1 = 4 ; 
sample_last_SEP_t1 = 4;

    for i= sample_first_SEP_t1:sample_last_SEP_t1
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('September t1');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t1',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.5 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end

% t-peaks vK diagram Oct t0
subplot(4,2,3);
FigTpeaksSEP.Color=([1 1 1]);
FigTpeaksSEP.Units='centimeters';
FigTpeaksSEP.Position=([5 3 50 25]);
dotsize = 50;
cut= 0;
sample_first_OCT_t0 = 9 ; 
sample_last_OCT_t0 = 11;

    for i= sample_first_OCT_t0:sample_last_OCT_t0
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            ylabel('H/C'),xlabel('O/C')
            title('October t0');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t0',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.07 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end


% t-peaks vK diagram Oct t1
subplot(4,2,4)
FigTpeaksOct.Color=([1 1 1]);
dotsize = 50;
cut= 0;
sample_first_OCT_t1 = 12 ; 
sample_last_OCT_t1 = 12;

    for i= sample_first_OCT_t1:sample_last_OCT_t1
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('October t1');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t1',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.5 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end


% t-peaks vK diagram December t0
subplot(4,2,5);
FigTpeaksSEP.Color=([1 1 1]);
FigTpeaksSEP.Units='centimeters';
FigTpeaksSEP.Position=([5 3 50 25]);
dotsize = 50;
cut= 0;
sample_first_DEC_t0 = 17 ; 
sample_last_DEC_t0 = 19;

    for i= sample_first_DEC_t0:sample_last_DEC_t0
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            ylabel('H/C'),xlabel('O/C')
            title('December t0');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t0',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.07 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end


% t-peaks vK diagram December t1
subplot(4,2,6)
FigTpeaksSEP.Color=([1 1 1]);
dotsize = 50;
cut= 0;
sample_first_DEC_t1 = 20 ; 
sample_last_DEC_t1 = 20;

    for i= sample_first_DEC_t1:sample_last_DEC_t1
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('December t1');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t1',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.5 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end

%February t-peaks
 
% t-peaks vK diagram Feb t0
subplot(4,2,7);
FigTpeaksSEP.Color=([1 1 1]);
FigTpeaksSEP.Units='centimeters';
FigTpeaksSEP.Position=([5 3 50 25]);
dotsize = 50;
cut= 0;
sample_first_FEB_t0 = 25 ; 
sample_last_FEB_t0 = 27;

    for i= sample_first_FEB_t0:sample_last_FEB_t0
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            ylabel('H/C'),xlabel('O/C')
            title('February t0');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t0',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.07 .5 .5 .5],'String','(a)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end


% t-peaks vK diagram Feb t1
subplot(4,2,8)
FigTpeaksFEB.Color=([1 1 1]);
dotsize = 50;
cut= 0;
sample_first_FEB_t1 = 28 ; 
sample_last_FEB_t1 = 28;

    for i= sample_first_FEB_t1:sample_last_FEB_t1
        plotposition = smpdat(i,1); %go to column  on smpdat and change to chronological order instead of just 1
        if (plotposition ~= 0) %change smpdat column 2 to chronoligal order in order to plot each sample i.e. 1.2.3...
            datplot=avg_intensities(:,i);
            a=datplot>cut;
            hold on
            scatter(formula(a&~Tpeaks_formula_filter,3),formula(a&~Tpeaks_formula_filter,2),dotsize,'k', 'filled') %Present only in data
            scatter(formula(a&Tpeaks_formula_filter,3),formula(a&Tpeaks_formula_filter,2),dotsize,'MarkerFaceColor','#FFC20A') %Present in paper and data
            axis([0 1 0 2.5])
            %ylabel('H/C')
            xlabel('O/C')
            title('February t1');
            sampling = char(sampcodes(i));
            month = sampling(1:3);
            legend('t1',strcat('t-Peaks = ',num2str(sum(a&Tpeaks_formula_filter))),'Location','bestoutside')
            set(gca,'FontSize', 19);
            %annotation('textbox',[.5 .5 .5 .5],'String','(b)','EdgeColor','none', 'FontSize', 24) %(x y width height)
        end 
    end
%% Tpeaks matches
%Tpeaks_formula_filter(2229,1) #Booleean
%avg_intensities(2229,49) Double
avg_intensities_logical = logical(avg_intensities);
avg_intensities_logical_Tpeaks = avg_intensities_logical(Tpeaks_formula_filter,:);
formulas_Tpeaks_matches = formula(Tpeaks_formula_filter,:);

number = size(avg_intensities,2);
tpeak_matches_i_and_j = zeros(number,number);
tpeak_matches_i_only = zeros(number,number);
tpeak_matches_j_only = zeros(number,number);

%make correlation map.
for i = 1:number 
    for j = 1:number
        tpeak_matches_i_and_j(i,j) = sum(avg_intensities_logical_Tpeaks(:,i) & avg_intensities_logical_Tpeaks(:,j));
        tpeak_matches_i_only(i,j) = sum(avg_intensities_logical_Tpeaks(:,i).* ~avg_intensities_logical_Tpeaks(:,j));
        tpeak_matches_j_only(i,j) = sum(avg_intensities_logical_Tpeaks(:,j).* ~avg_intensities_logical_Tpeaks(:,i));
    end
end

%How to get numbers after running script
% tpeak_matches_i_and_j(1,4) % Number Tpeaks present in both sample 1 and 4
% tpeak_matches_i_only(1,4) % Number of Tpeaks present in sample 1 and not 4
% tpeak_matches_i_only(4,1) % Number of Tpeaks present in sample 4 and not1
% tpeak_matches_i_and_j(1,1) % Number Tpeaks in sample 1
%tpeak_matches_i_and_j(2,2) % Number Tpeaks in sample 2
%tpeak_matches_i_and_j(3,3) % Number Tpeaks in sample 3

%% T-peaks ttest
%September
if true
%Student T-test (unpaired/independent)  h=1 rejects null hypothesis at 5% significance level (i.e. 95% confidence level) 
tpeaks_sept0 = [tpeak_matches_i_and_j(1,1), tpeak_matches_i_and_j(2,2), tpeak_matches_i_and_j(3,3)];
tpeaks_sept1 = [tpeak_matches_i_and_j(4,4), tpeak_matches_i_and_j(5,5), tpeak_matches_i_and_j(6,6)];
[h, p, ci, stats] = ttest2(tpeaks_sept0, tpeaks_sept1);
% Display results
disp(['t-value September: ', num2str(stats.tstat)]);
disp(['Degrees of Freedom September: ', num2str(stats.df)]);
disp(['p-value September: ', num2str(p)]);
 
%October
%Student T-test (unpaired/independent)  h=1 rejects null hypothesis at 5% significance level (i.e. 95% confidence level) 
tpeaks_octt0 = [tpeak_matches_i_and_j(9,9), tpeak_matches_i_and_j(10,10), tpeak_matches_i_and_j(11,11)];
tpeaks_octt1 = [tpeak_matches_i_and_j(12,12), tpeak_matches_i_and_j(13,13), tpeak_matches_i_and_j(14,14)];
[h, p, ci, stats] = ttest2(tpeaks_octt0, tpeaks_octt1);
% Display results
disp(['t-value October: ', num2str(stats.tstat)]);
disp(['Degrees of Freedom October: ', num2str(stats.df)]);
disp(['p-value October: ', num2str(p)]);

%December
%Student T-test (unpaired/independent)  h=1 rejects null hypothesis at 5% significance level (i.e. 95% confidence level) 
tpeaks_dect0 = [tpeak_matches_i_and_j(17,17),tpeak_matches_i_and_j(18,18), tpeak_matches_i_and_j(19,19)];
tpeaks_dect1 = [tpeak_matches_i_and_j(20,20), tpeak_matches_i_and_j(21,21), tpeak_matches_i_and_j(22,22)];
[h, p, ci, stats] = ttest2(tpeaks_dect0, tpeaks_dect1);
% Display results
disp(['t-value December: ', num2str(stats.tstat)]);
disp(['Degrees of Freedom December: ', num2str(stats.df)]);
disp(['p-value December: ', num2str(p)]);

%February
%Student T-test (unpaired/independent)  h=1 rejects null hypothesis at 5% significance level (i.e. 95% confidence level) 
tpeaks_febt0 = [tpeak_matches_i_and_j(25,25), tpeak_matches_i_and_j(26,26), tpeak_matches_i_and_j(27,27)];
tpeaks_febt1 = [tpeak_matches_i_and_j(28,28), tpeak_matches_i_and_j(29,29), tpeak_matches_i_and_j(30,30)];
[h, p, ci, stats] = ttest2(tpeaks_febt0, tpeaks_febt1);
% Display results
disp(['t-value February: ', num2str(stats.tstat)]);
disp(['Degrees of Freedom February: ', num2str(stats.df)]);
disp(['p-value February: ', num2str(p)]);
end 


%% Table DOM metrics
form_rnd=round(formula(:,8),4);
weights=avg_intensities(:,smp);
metricvalues_HC= formula(:,2);
for smp=1:nosamps
    nopeaks(smp,1)=sum(avg_intensities(:,smp)>0);
    tot_int=sample_intensities{smp}(smp_confident(:,smp),:);
    TAC(smp,1)=sum(tot_int(:));
    %effective_sample_size=sum(avg_intensities(:,smp))^2 / sum(avg_intensities(:,smp).^2);
    effective_sample_size=nopeaks(smp);
    OC_weighted_average=sum(avg_intensities(:,smp).*formula(:,3))/sum(avg_intensities(:,smp));
    OC_wa(smp,1)=OC_weighted_average;
    OC_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,3)-OC_weighted_average).^2)/sum(avg_intensities(:,smp));
    OC_wvar(smp,1)=OC_weighted_variance;
    OC_weighted_std_deviation = sqrt(OC_weighted_variance);
    OC_wstdev(smp,1) = OC_weighted_std_deviation ;
    OC_standard_error_mean=OC_weighted_std_deviation/sqrt(effective_sample_size);
    OC_SEM(smp,1) = OC_standard_error_mean;
    HC_weighted_average=sum(avg_intensities(:,smp).*formula(:,2))/sum(avg_intensities(:,smp));
    HC_wa(smp,1)=HC_weighted_average;
    HC_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,2)-HC_weighted_average).^2)/sum(avg_intensities(:,smp));
    HC_wvar(smp,1)=HC_weighted_variance;
    HC_weighted_std_deviation = sqrt(HC_weighted_variance);
    HC_wstdev(smp,1) = HC_weighted_std_deviation ;
    HC_standard_error_mean=HC_weighted_std_deviation/sqrt(effective_sample_size);
    HC_SEM(smp,1) = HC_standard_error_mean;
    MW_weighted_average=sum(avg_intensities(:,smp).*formula(:,1))/sum(avg_intensities(:,smp));
    MW_wa(smp,1)=MW_weighted_average;
    MW_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,1)-MW_weighted_average).^2)/sum(avg_intensities(:,smp));
    MW_wvar(smp,1)=MW_weighted_variance;
    MW_weighted_std_deviation = sqrt(MW_weighted_variance);
    MW_wstdev(smp,1) = MW_weighted_std_deviation ;
    MW_standard_error_mean=MW_weighted_std_deviation/sqrt(effective_sample_size);
    MW_SEM(smp,1) = MW_standard_error_mean;
    %
    cut = 0; 
    datplot=avg_intensities(:,smp); %plot avg intensities. loop can only use i 
    %
    c=datplot>cut&elements_used(:,4)>0&elements_used(:,5)==0; % requires N>0 S=0, CHON
    CHON(smp,1)=sum(c);
    CHON_relintensity(smp,1) = sum(datplot(c))./sum(datplot);
    
    b=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==0;% N=0 S=0,  CHO
    CHO(smp,1)=sum(b);
    CHO_relintensity(smp,1) = sum(datplot(b))./sum(datplot); %gives the decimal for average relative abundance of CHO in each sample (i.e. take sum of CHO intensities and divide by the total abundance of CHO for each sample)
    
    
    d=datplot>cut&elements_used(:,4)==0&elements_used(:,5)==1; %N=0 S=1, CHOS
    CHOS(smp,1)=sum(d);
    CHOS_relintensity(smp,1) = sum(datplot(d))./sum(datplot);
    DBE(smp,1)=sum(avg_intensities(:,smp).*formula(:,7))/sum(avg_intensities(:,smp));% DBE weighted average col 9 
    DBEminO(smp,1)= sum(avg_intensities(:,smp).*formula(:,8))/sum(avg_intensities(:,smp));% DBE-O weighted average col 10 
    DBEminO_wstdev(smp,1) = std(formula(:,8),avg_intensities(:,smp));
    AImod_weighted_average=sum(avg_intensities(:,smp).*formula(:,10))/sum(avg_intensities(:,smp));
    AImod_wa(smp,1)=AImod_weighted_average;
    AImod_weighted_variance=sum(avg_intensities(:,smp).*(formula(:,10)-AImod_weighted_average).^2)/sum(avg_intensities(:,smp));
    AImod_wvar(smp,1)=AImod_weighted_variance;
    AImod_weighted_std_deviation = sqrt(AImod_weighted_variance);
    AImod_wstdev(smp,1) = AImod_weighted_std_deviation ;
    AImod_standard_error_mean=AImod_weighted_std_deviation/sqrt(effective_sample_size);
    AImod_SEM(smp,1) = AImod_standard_error_mean;
    
    
end
clear  b c d cut 
table2=[nopeaks TAC OC_wa HC_wa MW_wa CHON CHO CHOS DBE DBEminO CHON_relintensity CHO_relintensity CHOS_relintensity  OC_wstdev HC_wstdev MW_wstdev AImod_wa AImod_wstdev DBEminO_wstdev HC_SEM OC_SEM AImod_SEM MW_SEM]; 
%% goodsmp badsmp
goodsmp=1:nosamps;
badsmp = smpdat(:,1)==0|table2(:,1)<700; 
goodsmp(badsmp)=[]; 
sampcodes_cut=sampcodes(goodsmp);

%% TAC TIC signals
for i=1:nosamps
    TICx(i,1)=sum(sample_data{i,1}(:,2));
    TICn(i,1)=sum(raw_data{i,1}(:,2));
end
figure();
subplot(2,1,1)
hold on
plot(1:nosamps,TAC)
plot(1:nosamps,TICx)
plot(1:nosamps,TICn)
legend('TAC','signal','TIC')
subplot(2,1,2)
hold on
plot(1:nosamps,TAC./TICx)
legend('TAC/signal')
%% Bray curtis dissimilarity 
dist_matrix=zeros(nosamps,nosamps); %creating a matrix template based on sample number
for dist_x=1:nosamps
    for dist_y=1:nosamps
        set1=avg_intensities(:,dist_x);
        set2=avg_intensities(:,dist_y);
        T1=set1/sum(set1); % T1 T2 are normalized here
        T2=set2/sum(set2);
        T3=abs(T1-T2); %
        T4=T1+T2;
        diss=sum(T3(:))/sum(T4(:))*100; %B-C Diss formula used
        dist_matrix(dist_x,dist_y)=diss;
    end
end
clear T1 T2 T3 T4 diss


dist_avg=mean(dist_matrix,1);
mean_samp_no=find(dist_avg==min(dist_avg));
dist_matrix_cut=dist_matrix(goodsmp,goodsmp);
BCD_U=triu(dist_matrix_cut);
b=find(BCD_U~=0);
median(BCD_U(b))
quantile(BCD_U(b),0.75)
quantile(BCD_U(b),0.25)
max(BCD_U(b))

[Y,eigvals] = cmdscale(dist_matrix_cut); %use Y eig values to make PCoA



%% dendrogram
%based on bray curtis dissimilarity values
sampcodes_cut=sampcodes(goodsmp);
FigDend = figure(); 
FigDend.Color=([1 1 1]);
FigDend.Units='centimeters';
FigDend.Position=([5 5 12 12]);
colormap winter
[rows,columns] = size(dist_matrix_cut);
v_mat = [];
    for i = 1:rows-1
                v_mat = [v_mat dist_matrix_cut(i,i+1:columns)];
    end
    Z=linkage(v_mat);
    Clu = cluster(Z,'Cutoff',15,'Criterion','distance');
[H,T,outperm]=dendrogram(Z,0,'ColorThreshold',15);
    axis([0.5 length(goodsmp)+0.5 0 50])
    xticklabels(sampcodes_cut(outperm))
    xtickangle(90)
    ticklabels = get(gca,'XTickLabel');
 lineColours = cell2mat(get(H,'Color'));  
 colourList = unique(lineColours, 'rows');

myColours = [230,186,0;
             127,127,127;
             10,35,140;
             176,0,0;
             158,182,72;
             79,129,189;
             209,99,9;
             4,122,146]/255;

%// Replace each colour (colour by colour). Start from 2 because the first colour are the "unclustered" black lines             
for colour = 2:size(colourList,1)
    %// Find which lines match this colour
    idx = ismember(lineColours, colourList(colour,:), 'rows');
    %// Replace the colour for those lines
    lineColours(idx, :) = repmat(myColours(colour-1,:),sum(idx),1);
end
%// Apply the new colours to the chart's line objects (line by line)
for linex = 1:size(H,1)
    set(H(linex), 'Color', lineColours(linex,:));
end 


%% formula presence 
% presence to absence of formulas (total 2229 formulas)
figure ('name','formulas:presence sorted')
mzthreshold = 0;
Ngoodsamples = length(goodsmp);
avg_intensities_logical = avg_intensities(:,goodsmp) > mzthreshold;
sumdata = sum(avg_intensities_logical,2); % 0 if present in no samples, Ngoodsamples if present in all samples, etc..
presentInAll = (sumdata == Ngoodsamples);
NpresentInAll = sum(presentInAll) ;%number of formulas present in all samples
plot(sumdata);
plot(sort(sumdata));


%% Figure 6 and S5 Histograms

% Histograms Filtered 
FigFilteredHistogram = figure('Name','Filtered Histograms');
FigFilteredHistogram.Position=([5 3 1500 800]);
hold on

fontsizeall=24;
figurelettersize=24;
linewidth=0.99;
histxlabel='H/C';
histylabelformulas='Number of formulas';
histylabelIntensities='Normalized intensities';
histylimmaxPeaks=690;
histylimmaxIntensities=500000;
histxlimmin=0.3;
histxlimmax=2.3;
numxticks=5;

colort0 = "#FFC20A";%yellow orange
colort1 ="#556B2F"; %dark green

%
t0 = {};
t0{1} = [1:3]; %September
t0{2} = [9:11]; %October
t0{3} = [17:19]; %December
t0{4} = [25:27]; %February

t1 = {};
t1{1} = [4]; %September
t1{2} = [12]; %October
t1{3} = [20]; %December
t1{4} = [28]; %February

month = {};
month{1} = 'September';
month{2} = 'October';
month{3} = 'December';
month{4} = 'February';

hist_HCratio = formula(:,2);
hist_data = avg_intensities;
hist_data_binary = hist_data ~= 0;

nbins = 8;
bins = linspace(histxlimmin,histxlimmax,nbins+1);
bin_centers = bins(1:end-1) + diff(bins) / 2;

%Peak number histograms (F)
for i=1:size(t0,2)
    subplot(2,size(t0,2),i)
    hist_data_t0 = sum(hist_data_binary(:,t0{i}),2)>0;
    hist_data_t1 = sum(hist_data_binary(:,t1{i}),2)>0;
    [counts_t0] = histcounts(hist_HCratio(hist_data_t0), bins);
    [counts_t1] = histcounts(hist_HCratio(hist_data_t1), bins);
    b=bar(bin_centers, [counts_t0;counts_t1]','grouped');
    b(1).FaceColor = colort0; 
    b(2).FaceColor = colort1; 
    ylim([0 histylimmaxPeaks]);
    xlim([histxlimmin histxlimmax]);
    num_ticks = numxticks; % Adjust this number for more or fewer x ticks
    xtickpositions=linspace(histxlimmin, histxlimmax, num_ticks);
    xticks(linspace(histxlimmin, histxlimmax, num_ticks));
    xticklabels(arrayfun(@(x) sprintf('%.1f', x), xtickpositions, 'UniformOutput', false));
    if i==1 
       ylabel(histylabelformulas);
    else
        set(gca,'YTickLabel',[])
    end  
    if i==4 
        legend('t0','t1')
    end
    set(gca,'FontSize', fontsizeall);
    title(month{i});
    set(gca, 'LineWidth', linewidth); % Thicker box around the plot area
    set(gca, 'Box', 'on'); 
end

% weighted intensities histograms (F)
for i=1:size(t0,2)
    bin_intensity_t0 = zeros(nbins,1);
    bin_intensity_t1 = zeros(nbins,1);
    weights_t0 = mean(avg_intensities(:,t0{i}),2);
    weights_t1 = mean(avg_intensities(:,t1{i}),2);
    for j=1:size(weights_t0,1)
        bin_n = find(bins < hist_HCratio(j),1,'last');
        if ~isempty(bin_n)
            bin_intensity_t0(bin_n) = bin_intensity_t0(bin_n) + weights_t0(j);
            bin_intensity_t1(bin_n) = bin_intensity_t1(bin_n) + weights_t1(j);
        end
    end
    subplot(2,size(t0,2),(i+4))
    hold on
    b=bar(bin_centers, [transpose(bin_intensity_t0);transpose(bin_intensity_t1)]','grouped');
    b(1).FaceColor = colort0; 
    b(2).FaceColor = colort1;
    ylim([0 histylimmaxIntensities]);
    xlim([histxlimmin histxlimmax]);
    num_ticks = numxticks; % Adjust this number for more or fewer x ticks
    xtickpositions=linspace(histxlimmin, histxlimmax, num_ticks);
    xticks(linspace(histxlimmin, histxlimmax, num_ticks));
    xticklabels(arrayfun(@(x) sprintf('%.1f', x), xtickpositions, 'UniformOutput', false));
    if i==1 
        ylabel(histylabelIntensities);
    else
        set(gca,'YTickLabel',[])
    end
    set(gca,'FontSize', fontsizeall); 
    set(gca, 'LineWidth', linewidth); % Thicker box around the plot area
    set(gca, 'Box', 'on'); 
end

% Histograms Unfiltered

FigUnfilteredHistogram = figure('Name','Unfiltered Histograms');
FigUnfilteredHistogram.Position=([5 3 1500 800]);
hold on

fontsizeall=24;
figurelettersize=24;
histxlabel='H/C';
histylabelformulas='Number of formulas';
histylabelIntensities='Normalized intensities';
histylimmaxPeaks=690;
histxlimmin=0.3;
histxlimmax=2.3;
numxticks=5;

colort0 = "#FFC20A";%yellow orange
colort1 ="#556B2F"; %dark green

t0 = {};
t0{1} = [1:3]; %September
t0{2} = [9:11]; %October
t0{3} = [17:19]; %December
t0{4} = [25:27]; %February

t1UF = {};
t1UF{1} = [7]; %September
t1UF{2} = [15]; %October
t1UF{3} = [23]; %December
t1UF{4} = [31]; %February

month = {};
month{1} = 'September';
month{2} = 'October';
month{3} = 'December';
month{4} = 'February';

nbins = 8;
bins = linspace(histxlimmin,histxlimmax,nbins+1);
bin_centers = bins(1:end-1) + diff(bins) / 2;

hist_HCratio = formula(:,2);
hist_data = avg_intensities;
hist_data_binary = hist_data ~= 0;

%Peak number histograms (UF)
for i=1:size(t0,2)
    subplot(2,size(t0,2),i)
    hist_data_t0 = sum(hist_data_binary(:,t0{i}),2)>0;
    hist_data_t1_UF = sum(hist_data_binary(:,t1UF{i}),2)>0;
    [counts_t0] = histcounts(hist_HCratio(hist_data_t0), bins);
    [counts_t1_UF] = histcounts(hist_HCratio(hist_data_t1_UF), bins);
    b=bar(bin_centers, [counts_t0;counts_t1_UF]','grouped');
    b(1).FaceColor = colort0; 
    b(2).FaceColor = colort1;   
    ylim([0 histylimmaxPeaks]);
    xlim([histxlimmin histxlimmax]);
    num_ticks = numxticks; % Adjust this number for more or fewer x ticks
    xtickpositions=linspace(histxlimmin, histxlimmax, num_ticks);
    xticks(linspace(histxlimmin, histxlimmax, num_ticks));
    xticklabels(arrayfun(@(x) sprintf('%.1f', x), xtickpositions, 'UniformOutput', false));
    if i==1 
       ylabel(histylabelformulas);
    else
        set(gca,'YTickLabel',[])
    end
    if i==4 
        legend('t0','t1UF')
    end
    set(gca, 'LineWidth', linewidth); % Thicker box around the plot area
    set(gca, 'Box', 'on');
    set(gca,'FontSize', fontsizeall);
    title(month{i});
end

% weighted intensities histograms (UF)

for i=1:size(t0,2)
    bin_intensity_t0 = zeros(nbins,1);
    bin_intensity_t1UF = zeros(nbins,1);
    weights_t0 = mean(avg_intensities(:,t0{i}),2);
    weights_t1UF = mean(avg_intensities(:,t1UF{i}),2);
    for j=1:size(weights_t0,1)
        bin_n = find(bins < hist_HCratio(j),1,'last');
        if ~isempty(bin_n)
            bin_intensity_t0(bin_n) = bin_intensity_t0(bin_n) + weights_t0(j);
            bin_intensity_t1UF(bin_n) = bin_intensity_t1UF(bin_n) + weights_t1UF(j);
        end
    end
    subplot(2,size(t0,2),(i+4))
    hold on
    b=bar(bin_centers, [transpose(bin_intensity_t0);transpose(bin_intensity_t1UF)]','grouped');
    b(1).FaceColor = colort0; 
    b(2).FaceColor = colort1;
    ylim([0 histylimmaxIntensities]);
    xlim([histxlimmin histxlimmax]);
    num_ticks = numxticks; % Adjust this number for more or fewer x ticks
    xtickpositions=linspace(histxlimmin, histxlimmax, num_ticks);
    xticks(linspace(histxlimmin, histxlimmax, num_ticks));
    xticklabels(arrayfun(@(x) sprintf('%.1f', x), xtickpositions, 'UniformOutput', false));
    if i==1 
        ylabel(histylabelIntensities);
    else
        set(gca,'YTickLabel',[])
    end
    set(gca,'FontSize', fontsizeall);
    set(gca, 'LineWidth', linewidth); % Thicker box around the plot area
    set(gca, 'Box', 'on');
end
%% Figure S9 SPE-DOC Recovery   
samplename = categorical (sampcodes_cut);
samplename = reordercats(samplename,string(samplename));

if true
    groupingvariable1 = smpdat(goodsmp,2); % column 2 is numbered for months and incubation need to remove 0s for goodsmp time
    %Boxplot for SPE-DOC % recovery
    FigBoxplotSPEDOC = figure('Name','SPE-DOC % recovery');
    boxplot(smpdat(goodsmp,3), groupingvariable1,'Labels',{'Sep t0','Sep t1','Oct t0','Oct t1','Dec t0','Dec t1','Feb t0','Feb t1'});
    ylabel('SPE-DOC % recovery');
end
%% assigned data
form_head={'C','H','O','N','S','P','Se','Na','Cl','13C','m/z','H/C','O/C','KM','NKM','KMD','DBE','DBE-O','13Cratio','AImod'};
for smp=1:nosamps
    headx{1,smp}=strcat(sampcodes{smp});
end
    heads=[form_head headx];
    assign_dat=[elements_used formula avg_intensities];
    csvwrite_with_headers('assigned_data.csv',assign_dat,heads)