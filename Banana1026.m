%% This script calculates iteratively erosion rates and burial times from 10Be and 26Al concentrations scaled to their catchment average production rates (samples.csv).
%% The distribution of spallation and muon relative production rates can be defined in the file procuctions.csv. They can be calculated using this approach: https://angelrodes.wordpress.com/2021/12/15/average-cosmogenic-production-rate-calculator/, or just with the latitude-elevation approximations defined in Rodés, Á. (2021) "The NUNAtak Ice Thinning (NUNAIT) Calculator for Cosmonuclide Elevation Profiles" Geosciences 11, no. 9: 362. doi:10.3390/geosciences11090362 
%% The graphical outputs are a banana plot and the results in the time-erosion rate spaces.
%% In the command window, results are shown as ranges derived from the concentration uncertainties.
%% Angel Rodes, 2022
%% www.angelrodes.com

clear
close all hidden
clc
%% Input samples
selectedfile = 'samples.csv';
fid = fopen(selectedfile);
mydata = textscan(fid, '%s %f %f %f %f ',...
    'HeaderLines', 1,'Delimiter',',');
fclose(fid);
samplename=mydata{1}';
C10=mydata{2}';
dC10=mydata{3}';
C26=mydata{4}';
dC26=mydata{5}';

%% Input productions
selectedfile = 'productions.csv';
fid = fopen(selectedfile);
mydata = textscan(fid, ' %f %f %f ',...
    'HeaderLines', 1,'Delimiter',',');
fclose(fid);
P10=mydata{1};
P26=mydata{2};
L=mydata{3};

if max([numel(P10),numel(P26),numel(L)])>4
    warning('Too many production rates. (max=4)')
end

if range([numel(P10),numel(P26),numel(L)])>0
    warning('Production rate lengths should be consistent.')
end

% make production lengths==4
for n=1:4
    if numel(L)<n
        P10(n)=0;
        P26(n)=0;
        L(n)=0;
    end
end


% Define half lives
l10=log(2)/1.387e6;
l26=log(2)/0.705e6;

% Define density
density=2.65;


for sample=1:numel(C10)

    %% Make e-t space
    % for plotting purposes (avoid non-zero concentrations), the sample has been considered to be
    % considered to be transported in one year, accumulating P*1 in the basin
    % t (a) and e (m/Ma)
    N10=@(er,t) (P10(1)+P10(2)+P10(3)+P10(4))+...
        exp(-l10.*t).*(...
        P10(1)./(l10+density.*er./1e4./L(1)) +...
        P10(2)./(l10+density.*er./1e4./L(2)) +...
        P10(3)./(l10+density.*er./1e4./L(3)) +...
        P10(4)./(l10+density.*er./1e4./L(4)) ...
            );
    N26=@(er,t) (P26(1)+P26(2)+P26(3)+P26(4))+...
        exp(-l26.*t).*(...
        P26(1)./(l26+density.*er./1e4./L(1)) +...
        P26(2)./(l26+density.*er./1e4./L(2)) +...
        P26(3)./(l26+density.*er./1e4./L(3)) +...
        P26(4)./(l26+density.*er./1e4./L(4)) ...
            );
    
    %% find solutions
    chisq=@(er,t) ((((N10(er,t)-C10(sample))./dC10(sample)).^2+((N26(er,t)-C26(sample))./dC26(sample)).^2));
    plotpoints=200;
    tref=[linspace(0,3E6,round(plotpoints/2)) linspace(3E6,2E7,round(plotpoints/2))];
    eref=logspace(-3,6,plotpoints);
    convergence=0;
    minsolutions=30;
    count=0;
    while convergence<1
        count=count+1;
        [tmatrix,ematrix] = meshgrid(tref,eref);
        chimatrix=chisq(ematrix,tmatrix);
        sel=find(chimatrix==min(chimatrix(:)));
        sel=sel(1);
        chimatrixbroken=chimatrix;
        chimatrixbroken(tmatrix==tmatrix(sel)|ematrix==ematrix(sel))=max(chimatrix(:)); % remove minimum rows and columns
        
        if min(chimatrixbroken(:))<1
            maxchi=min(chimatrixbroken(:))*2+2;
            selt=(min(chimatrix,[],1)<=maxchi);
            tref=tref(selt);
            tref=interp1(tref,linspace(1,length(tref),plotpoints));
            tref(1)=max(0,tref(1)-(tref(2)-tref(1)));
            tref(end)=tref(end)+(tref(end)-tref(end-1));
            sele=(min(chimatrix,[],2)<=maxchi);
            eref=eref(sele');
            eref=interp1(eref,linspace(1,length(eref),plotpoints));
            eref(1)=max(0,eref(1)-(eref(2)-eref(1)));
            eref(end)=eref(end)+(eref(end)-eref(end-1));
        else
            selt=(min(chimatrix,[],1)<=max(2,prctile(min(chimatrix,[],1),50)));
            tref=tref(selt);
            tref=interp1(tref,linspace(1,length(tref),plotpoints));
            tref(1)=max(0,tref(1)-(tref(2)-tref(1)));
            tref(end)=tref(end)*(tref(end)/tref(end-1));
            sele=(min(chimatrix,[],2)<=max(2,prctile(min(chimatrix,[],2),50)));
            eref=eref(sele');
            eref=interp1(eref,linspace(1,length(eref),plotpoints));
            eref(1)=max(0,eref(1)-(eref(2)-eref(1)));
            eref(end)=eref(end)*(eref(end)/eref(end-1));
        end
        
        
        if sum(sum(chimatrix<1))>minsolutions
            convergence=convergence+0.2;
            [tmatrix,ematrix] = meshgrid(tref,eref);
            chimatrix=chisq(ematrix,tmatrix);
            chis=sort(unique(chimatrix(:)));
            maxchi=min(chis(chis>min(chimatrix(:))+1));
            sel=(chimatrix<=maxchi);
            sel2=find(chimatrix==min(chimatrix(:)));
            sel2=sel2(1);
            solutione=[unique(min(ematrix(sel))),unique(max(ematrix(sel))),ematrix(sel2)];
            solutiont=[unique(min(tmatrix(sel))),unique(max(tmatrix(sel))),tmatrix(sel2)];
            eprecis=round(log10(min(diff(unique(eref)))));
            tprecis=round(log10(min(diff(unique(tref)))));
            solutione=round(solutione,-eprecis);
            solutiont=round(solutiont,-tprecis);
            selmatrix=find(sel);
        end
        if count>1000
            convergence=1;
            solutione=[-1,-1,-1];
            solutiont=[-1,-1,-1];
        end
    end
    if sample==1
        disp(['Sample' '	' 'iterations'])
        disp(['______' '	' '__________'])
    end
        disp([samplename{sample} '	' num2str(count)])
    
    


   
    
    %% plot banana
    set(0, 'DefaultLineLineWidth', 2);
    
    if sample==1
        figure(1)
%         set(h,'Visible', 'off');
        %     figure(1)
        %     subplot(round(numel(C10)/2),2,sample)
        hold on
        plotpoints=1000;
        
%         % solutions
%         if solutione(1)>-1
%             plot(N10(ematrix(selmatrix),tmatrix(selmatrix)),N26(ematrix(selmatrix),tmatrix(selmatrix))./N10(ematrix(selmatrix),tmatrix(selmatrix)),...
%                 '.','Color',[0.5 0.5 0.5])
%         end
        
        % erosions
        ttemp=[0, logspace(0,8,plotpoints), Inf];
        for etemp=[0 1 10 100 1000 Inf]
            N10temp=N10(etemp,ttemp);
            N26temp=N26(etemp,ttemp);
            plot(N10temp,N26temp./N10temp,'-','Color',[0.5 0.5 0.5])
            sel=1;
            text(N10temp(sel),N26temp(sel)./N10temp(sel),[num2str(etemp) ' m/Ma'],...
                'Color',[0.5 0.5 0.5],'HorizontalAlignment','left','VerticalAlignment','bottom', 'Clipping', 'on')
        end
        
        % ages long
        etemp=[0,logspace(-2,7,plotpoints)];
        for ttemp=[0 0.5 1 2 3 4 5]*1e6
            N10temp=N10(etemp,ttemp);
            N26temp=N26(etemp,ttemp);
            plot(N10temp,N26temp./N10temp,'-','Color',[0.5 0.5 0.5])
            sel=1;
            text(N10temp(sel),N26temp(sel)./N10temp(sel),[num2str(ttemp/1e6) ' Ma'],...
                'Color',[0.5 0.5 0.5],'HorizontalAlignment','left','VerticalAlignment','top', 'Clipping', 'on')
        end
        
%         % ages short
%         etemp=[0,logspace(-2,7,plotpoints)];
%         for ttemp=[0 1 10 50 100 200 300 400 Inf]*1e3
%             N10temp=N10(etemp,ttemp);
%             N26temp=N26(etemp,ttemp);
%             plot(N10temp,N26temp./N10temp,':b')
%             sel=length(N10temp);
%             text(N10temp(sel),N26temp(sel)./N10temp(sel),[num2str(ttemp/1e3) ' ka'],...
%                 'Color','b','HorizontalAlignment','right','VerticalAlignment','top', 'Clipping', 'on')
%         end
    end
   
    % plot sample
    figure(1)
    hold on
    plotpoints=100;
    theta=linspace(0,2*pi,plotpoints);
    C10theta=C10(sample)+sin(theta)*dC10(sample);
    C26theta=C26(sample)+cos(theta)*dC26(sample);
    xaxissample=C10theta;
    yaxissample=C26theta./C10theta;
    plot(xaxissample,yaxissample,'Color','k','LineWidth',1)
    plot(C10(sample),C26(sample)/C10(sample),'*k')
    
    xlabel('[^{10}Be]')
    ylabel('[^{26}Al] / [^{10}Be]')
    
    
    text(C10(sample),C26(sample)/C10(sample)+dC26(sample)/C10(sample),samplename{sample},'HorizontalAlignment','center')
    
    set(gca,'xscale','log')
    box on
    xlim([min(C10-dC10)/2 N10(0,0)*2])
    ylim([0 max(C26./C10*1.2)])
    
    
    %% store results and calculate corrections
    if sample==1
        strsolutionsheader=['Sample_name' '	' ...
            'tmin(Ma)' '	' 'tmax(Ma)' '	' 'tbest(Ma)' '	' 'emin(m/Ma)' '	' 'emax(m/Ma)'  '	' 'ebest(m/Ma)'  '	' 'N'];
        %         strsolutionsheader=['Sample_name' '	' 'P10sp(basin)' '	' 'P10mu(basin)' '	' 'P10sp(post-burial)' '	' 'P10mu(post-burial)' '	'...
        %             'P26sp(basin)' '	' 'P26mu(basin)' '	' 'P26sp(post-burial)' '	' 'P26mu(post-burial)' '	'...
        %             'tmin(Ma)NP' '	' 'tmax(Ma)NP' '	' 'tbest(Ma)NP' '	' 'emin(m/Ma)NP' '	' 'emax(m/Ma)NP'  '	' 'ebest2(m/Ma)NP' '	' 'max % of muon produced' '	' '(NP=No post burial production considered)'];
        
    end
    
    
    strsolutions{sample}=[samplename{sample} '	'...
        num2str(solutiont(1)/1e6) '	' num2str(solutiont(2)/1e6) '	' num2str(solutiont(3)/1e6) '	' num2str(solutione(1)) '	' num2str(solutione(2)) '	' num2str(solutione(3)) '	' num2str(numel(selmatrix))];
    
        %% plot e,t resutls
        figure(2)
        hold on
        %     plot(solutiont(3)+1,solutione(3),'xb')
        %     plot([solutiont(1),solutiont(2)]+1,[solutione(3),solutione(3)],'-b')
        %     plot([solutiont(3),solutiont(3)]+1,[solutione(1),solutione(2)],'-b')
        contour(tmatrix+1,ematrix,chimatrix,[0,maxchi],'Color','b')
    sel=find(chimatrix<=maxchi);
    %     plot(tmatrix(sel),ematrix(sel),'.b')
        text(solutiont(3)+1,solutione(3),...
            samplename{sample},'Color','k','HorizontalAlignment','center','VerticalAlignment','middle', 'Clipping', 'on')
        set(gca,'yscale','log')
        % add 1 year to all ages if you use a log x scale
        % set(gca,'xscale','log')
        box on
        xlabel('Burial age')
        ylabel('Erosion rate (m/Ma)')
        title('Burial model results')
    
end

%% display results
clc
disp(strsolutionsheader)
for sample=1:numel(C10)
    disp(strsolutions{sample})
end
disp(' ')
disp('Angel Rodes, 2022')
disp('www.angelrodes.com')