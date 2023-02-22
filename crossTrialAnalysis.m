classdef crossTrialAnalysis < handle
    % Filename:	crossTrialAnalysis
    % Author:  James Hermus
    % Date:		1 April 2020
    % Description:	Object oriented program which calls the trialAnalysis
    % object program to computed pooled quanties across all 22 trials.
    
    properties
        subjDex
        speedDex
        trialDex
        stiffDex
        test % Structure from testAnalysis
        outputt_lnr_PC
        outputt_cv_v
        colorr
        
    end
    
    methods
        function this = crossTrialAnalysis(subjDex, speedDex, trialDex, stiffDex)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            if(nargin==0)
                this.subjDex = 1:3;
                this.speedDex = 1:2;
                this.trialDex = 1:3;
                this.stiffDex = 1;
            else
                this.subjDex = subjDex;
                this.speedDex = speedDex;
                this.trialDex = trialDex;
                this.stiffDex = stiffDex;
            end
            
            this.colorr = {'#4DBEEE','#A2142F'};
            
            this.get_testStruct();
            
        end
        
        function [] = get_testStruct(this)
            
            speed_StiffScale = [0.5,0.5,1.5,0.5,0.5,1.5]; % Scale for stiff conditions
            speed_DampingScale = [0.05,0.05,0.1,0.05,0.05,0.1]; % Scale for stiff conditions
            speed_Z_gain = [speed_StiffScale; speed_DampingScale]';
            
            stiff_scaleVec = [0.5000,    0.5000;...
                0.5000,    1.0000;...
                0.5000,    1.5000;...
                1.0000,    0.5000;...
                1.0000,    1.0000;...
                1.0000,    1.5000;...
                1.5000,    0.5000;...
                1.5000,    1.0000;...
                1.5000,    1.5000];
            
            stiff = this.stiffDex;
            clear test
            for subj = this.subjDex
                for speed = this.speedDex
                    Z_gain_input = stiff_scaleVec(stiff,:).*speed_Z_gain(speed,:);
                    for trial = this.trialDex
                        test{subj,speed,trial} = trialAnalysis(subj, speed, trial, Z_gain_input);
                    end
                end
                disp(['Subject: ',int2str(subj),', Stiffness Condition: ',int2str(stiff)]);
            end
                        
            this.test = test;
                        
        end
                
        function [v1,outputt_lnr_PC] = get_lnEig(this,lnr0,lnr45,anovaPlot)
            
            lnr0_meanTrial = mean(lnr0,1);
            lnr45_meanTrial = mean(lnr45,1);

            X = [lnr0_meanTrial(:),lnr45_meanTrial(:)];
            [V,D] = eig(X'*X);
            
            D = diag(D);
            
            [tmp,dex] = sort(D,'descend');
            
            D = D(dex);
            V = V(:,dex);
            v1 = V(:,1);
            v2 = V(:,2);
            
            % Sanity check
%             figure; plot([0,v1(1)],[0,v1(2)]); axis equal;

            % Compute PC1
            [n_trial, n_speed, n_dir, n_subj] = size(lnr0);
            for subj = 1:n_subj
                for speed = 1:n_speed
                    for dir = 1:n_dir
                        for trial = 1:n_trial
                            lnr_PC(trial,speed,dir,subj) = v1'*[lnr0(trial,speed,dir,subj);lnr45(trial,speed,dir,subj)];
                            if(dir == 1)
                                signFlip = 1;
                            elseif(dir == 2)
                                signFlip = -1;
                            end
                            lnr_PC_signFlip(trial,speed,dir,subj) = v1'*[signFlip.*lnr0(trial,speed,dir,subj);signFlip.*lnr45(trial,speed,dir,subj)];
                        end
                    end
                end
            end
            
            %% Run ANOVA on PC coordinates
            [outputt_lnr_PC] = this.threeWayANOVAJH(lnr_PC);
            [outputt_lnr_PC_signFlip] = this.threeWayANOVAJH(lnr_PC_signFlip);
 
%             lnr_PC_subjAve = squeeze(nanmean(lnr_PC,1));
%             [H,P,CI] = ttest(lnr_PC_subjAve(1,1,:),lnr_PC_subjAve(1,2,:)) % Slow
%             [H,P,CI] = ttest(lnr_PC_subjAve(2,1,:),lnr_PC_subjAve(2,2,:)) % Medium
%             [H,P,CI] = ttest(lnr_PC_subjAve(3,1,:),lnr_PC_subjAve(3,2,:)) % Fast
            
            if(anovaPlot)
                % Plot ln0
                figure; hold on;
                errorbar([0.075,0.5,2.0], outputt_lnr_PC.ave_CW, outputt_lnr_PC.dev_CW,'color',this.colorr{1},'linewidth',2.5);
                errorbar([0.075,0.5,2.0], outputt_lnr_PC.ave_CCW,outputt_lnr_PC.dev_CCW,'color',this.colorr{2},'linewidth',2.5);
%                 plot([0.075,0.5,2.0], outputt_lnr_PC.ave_CW,'color',this.colorr{1},'markersize',30);
%                 plot([0.075,0.5,2.0], outputt_lnr_PC.ave_CCW,'color',this.colorr{2},'markersize',30);
                            ylim([-0.7 0.7]); yticks(gca,[-0.6:0.3:0.6]);
                            xlim([0 2.075]); xticks(gca,[0,0.5,2]);

                for i = 1:10
                    plot([0.075],mean(lnr_PC(:,1,1,i)),'o','color',this.colorr{1});
                    plot([0.5],  mean(lnr_PC(:,2,1,i)),'o','color',this.colorr{1});
                    plot([2.0],  mean(lnr_PC(:,3,1,i)),'o','color',this.colorr{1});
                    plot([0.075],mean(lnr_PC(:,1,2,i)),'o','color',this.colorr{2});
                    plot([0.5],  mean(lnr_PC(:,2,2,i)),'o','color',this.colorr{2});
                    plot([2.0],  mean(lnr_PC(:,3,2,i)),'o','color',this.colorr{2});
                end
                ylabel('ln(r)_{PC}','fontsize',18); %ylim([0 0.4]);
                xlabel('Speed (rev/s)','fontsize',18);
                legend('CW','CCW','Location','North'); 
                set(gca,'FontSize',20,'linewidth',2.5);
                hold off;

                outputt_lnr_PC.P
                
                % Plot ln0 SignFlip
                figure; hold on;
                errorbar([0.075,0.5,2.0], outputt_lnr_PC_signFlip.ave_CW, outputt_lnr_PC_signFlip.dev_CW,'color',this.colorr{1},'linewidth',2.5);
                errorbar([0.075,0.5,2.0], outputt_lnr_PC_signFlip.ave_CCW,outputt_lnr_PC_signFlip.dev_CCW,'color',this.colorr{2},'linewidth',2.5);
                            ylim([-0.7 0.7]); yticks(gca,[-0.6:0.3:0.6]);
                            xlim([0 2.075]); xticks(gca,[0,0.5,2]);
                            
                for i = 1:10
                    plot([0.075],mean(lnr_PC_signFlip(:,1,1,i)),'o','color',this.colorr{1});
                    plot([0.5],  mean(lnr_PC_signFlip(:,2,1,i)),'o','color',this.colorr{1});
                    plot([2.0],  mean(lnr_PC_signFlip(:,3,1,i)),'o','color',this.colorr{1});
                    plot([0.075],mean(lnr_PC_signFlip(:,1,2,i)),'o','color',this.colorr{2});
                    plot([0.5],  mean(lnr_PC_signFlip(:,2,2,i)),'o','color',this.colorr{2});
                    plot([2.0],  mean(lnr_PC_signFlip(:,3,2,i)),'o','color',this.colorr{2});
                end

                ylabel('ln(r)_{PC} (CCW reflected)','fontsize',18); %ylim([0 0.4]);
                xlabel('Speed (rev/s)','fontsize',18);
                legend('CW','CCW','Location','North'); 
                set(gca,'FontSize',20,'linewidth',2.5);
                hold off;
            end
            
        end
        
        function [] = get_lnr0lnr45(this,lnr0,lnr45,PC1,subjectPlot)
            
              % Look at individual subjects
              if(subjectPlot)
                  for subj = this.subjDex
                      figure;
                      plot(squeeze(lnr0(:,:,1,subj)),squeeze(lnr45(:,:,1,subj)),'.','markersize',15); hold on;
                      plot(squeeze(lnr0(:,:,2,subj)),squeeze(lnr45(:,:,2,subj)),'.','markersize',15);
                      xlim([-1 1]); ylim([-1 1]); axis equal;
                      set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
                      
                      lr0_d1_mean = mean(lnr0(:,:,1,subj),1);
                      lr45_d1_mean = mean(lnr45(:,:,1,subj),1);
                      lr0_d1_std = std(lnr0(:,:,1,subj),0,1);
                      lr45_d1_std = std(lnr45(:,:,1,subj),0,1);
                      
                      lr0_d2_mean = mean(lnr0(:,:,2,subj),1);
                      lr45_d2_mean = mean(lnr45(:,:,2,subj),1);
                      lr0_d2_std = std(lnr0(:,:,2,subj),0,1);
                      lr45_d2_std = std(lnr45(:,:,2,subj),0,1);
                      
                      errorbar(lr0_d1_mean,lr45_d1_mean,...
                          lr45_d1_std,lr45_d1_std,...
                          lr0_d1_std,lr0_d1_std,'o','linewidth',2.5);
                      
                      errorbar(lr0_d2_mean,lr45_d2_mean,...
                          lr45_d2_std,lr45_d2_std,...
                          lr0_d2_std,lr0_d2_std,'o','linewidth',2.5);
                  end
              end
            
            for subj = this.subjDex
                lnr0_d1_trailMean(:,subj) = mean(lnr0(:,:,1,subj),1);
                lnr45_d1_trailMean(:,subj) = mean(lnr45(:,:,1,subj),1);
                lnr0_d2_trailMean(:,subj) = mean(lnr0(:,:,2,subj),1);
                lnr45_d2_trailMean(:,subj) = mean(lnr45(:,:,2,subj),1);
            end
            
            figure;
%             plot(lnr0_d1_trailMean',lnr45_d1_trailMean','.','markersize',25); hold on;
%             plot(lnr0_d2_trailMean',lnr45_d2_trailMean','.','markersize',25);
            markersizee = 50;
            ax1 = scatter(lnr0_d1_trailMean(1,:),lnr45_d1_trailMean(1,:),markersizee,'MarkerFaceColor',this.colorr{1},'MarkerEdgeColor',this.colorr{1},'MarkerFaceAlpha',0.3); hold on;
            ax2 = scatter(lnr0_d1_trailMean(2,:),lnr45_d1_trailMean(2,:),markersizee,'MarkerFaceColor',this.colorr{1},'MarkerEdgeColor',this.colorr{1},'MarkerFaceAlpha',0.6); hold on;
            ax3 = scatter(lnr0_d1_trailMean(3,:),lnr45_d1_trailMean(3,:),markersizee,'MarkerFaceColor',this.colorr{1},'MarkerEdgeColor',this.colorr{1},'MarkerFaceAlpha',1.0); hold on;

            ax4 = scatter(lnr0_d2_trailMean(1,:),lnr45_d2_trailMean(1,:),markersizee,'MarkerFaceColor',this.colorr{2},'MarkerEdgeColor',this.colorr{2},'MarkerFaceAlpha',0.3); hold on;
            ax5 = scatter(lnr0_d2_trailMean(2,:),lnr45_d2_trailMean(2,:),markersizee,'MarkerFaceColor',this.colorr{2},'MarkerEdgeColor',this.colorr{2},'MarkerFaceAlpha',0.6); hold on;
            ax6 = scatter(lnr0_d2_trailMean(3,:),lnr45_d2_trailMean(3,:),markersizee,'MarkerFaceColor',this.colorr{2},'MarkerEdgeColor',this.colorr{2},'MarkerFaceAlpha',1.0); 

            xlim([-0.75 0.75]); ylim([-0.75 0.75]); 
            set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin');
            
            lr0_d1_mean = mean(lnr0_d1_trailMean,2);
            lr45_d1_mean = mean(lnr45_d1_trailMean,2);
            lr0_d1_std = std(lnr0_d1_trailMean,0,2);
            lr45_d1_std = std(lnr45_d1_trailMean,0,2);
            
            lr0_d2_mean = mean(lnr0_d2_trailMean,2);
            lr45_d2_mean = mean(lnr45_d2_trailMean,2);
            lr0_d2_std = std(lnr0_d2_trailMean,0,2);
            lr45_d2_std = std(lnr45_d2_trailMean,0,2);
            
            e1 = errorbar(lr0_d1_mean(1),lr45_d1_mean(1),lr45_d1_std(1),lr45_d1_std(1),lr0_d1_std(1),lr0_d1_std(1),'color',this.colorr{1},'linewidth',2.5); 
            set([e1.Bar], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*0.3]);

            e2 = errorbar(lr0_d1_mean(2),lr45_d1_mean(2),lr45_d1_std(2),lr45_d1_std(2),lr0_d1_std(2),lr0_d1_std(2),'color',this.colorr{1},'linewidth',2.5); 
            set([e2.Bar], 'ColorType', 'truecoloralpha', 'ColorData', [e2.Line.ColorData(1:3); 255*0.6])
            
            errorbar(lr0_d1_mean(3),lr45_d1_mean(3),lr45_d1_std(3),lr45_d1_std(3),lr0_d1_std(3),lr0_d1_std(3),'color',this.colorr{1},'linewidth',2.5);

            e4 = errorbar(lr0_d2_mean(1),lr45_d2_mean(1),lr45_d2_std(1),lr45_d2_std(1),lr0_d2_std(1),lr0_d2_std(1),'color',this.colorr{2},'linewidth',2.5);
            set([e4.Bar, e4.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e4.Line.ColorData(1:3); 255*0.3])

            e5 = errorbar(lr0_d2_mean(2),lr45_d2_mean(2),lr45_d2_std(2),lr45_d2_std(2),lr0_d2_std(2),lr0_d2_std(2),'color',this.colorr{2},'linewidth',2.5);
            set([e5.Bar, e5.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e5.Line.ColorData(1:3); 255*0.6])

            errorbar(lr0_d2_mean(3),lr45_d2_mean(3),lr45_d2_std(3),lr45_d2_std(3),lr0_d2_std(3),lr0_d2_std(3),'color',this.colorr{2},'linewidth',2.5);

            
            ax_PC = plot(PC1(1)*[-1,1],PC1(2)*[-1,1],'--k','linewidth',2.5);
            set(gca, 'linewidth',2.5,'fontsize',16,'Box','off');
            xlabel('ln(r_0)'); ylabel('ln(r_{45})'); axis equal; xticks([-0.6 -0.3 0.3 0.6]); yticks([-0.6 -0.3 0.3 0.6]);
%             legend([ax1,ax2,ax3,ax4,ax5,ax6,ax_PC],...
%                 {'CW, slow','CW, Medium','CW, Fast','CCW, slow','CCW, Medium','CCW, Fast','PC_1'},...
%                 'location','southeast');

            
            %% Run ANOVA on ln0 and ln45
%                 
%             [outputt_ln0] = this.threeWayANOVAJH(lnr0);
%             [outputt_ln45] = this.threeWayANOVAJH(lnr45);
% 
%                         
%             % Plot ln0
%             figure; hold on;
%             errorbar([0.075,0.5,2.0], outputt_ln0.ave_CW, outputt_ln0.dev_CW,this.colorr(1),'linewidth',2.5);
%             errorbar([0.075,0.5,2.0], outputt_ln0.ave_CCW,outputt_ln0.dev_CCW,this.colorr(2),'linewidth',2.5);
%             plot([0.075,0.5,2.0], outputt_ln0.ave_CW,'.b','markersize',30);
%             plot([0.075,0.5,2.0], outputt_ln0.ave_CCW,'.r','markersize',30);
%             %             ylim([100 225]); yticks(gca,[135:45:225]);
%             %             xlim([0 2.25]); xticks(gca,[0:0.5:2]);
%             
% %             for i = 1:10
% %                 plot([0.5],lnr0(:,1,1,i),'.b');
% %                 plot([2.0],lnr0(:,2,1,i),'.b');
% %                 plot([0.5],lnr0(:,1,2,i),'.r');
% %                 plot([2.0],lnr0(:,2,2,i),'.r');
% %             end
%             ylabel('ln(r_{0})','fontsize',18); %ylim([0 0.4]);
%             xlabel('Speed (rev/s)','fontsize',18);
%             legend('CW','CCW','Location','North'); grid on;
%             set(gca,'FontSize',20);
%             hold off;
%             
%             figure;
%             hold on;
%             errorbar([0.075,0.5,2.0], outputt_ln45.ave_CW, outputt_ln45.dev_CW,this.colorr(1),'linewidth',2.5);
%             errorbar([0.075,0.5,2.0], outputt_ln45.ave_CCW,outputt_ln45.dev_CCW,this.colorr(2),'linewidth',2.5);
%             plot([0.075,0.5,2.0], outputt_ln45.ave_CW,'.b','markersize',30);
%             plot([0.075,0.5,2.0], outputt_ln45.ave_CCW,'.r','markersize',30);
%             %             ylim([100 225]); yticks(gca,[135:45:225]);
%             %             xlim([0 2.25]); xticks(gca,[0:0.5:2]);
%             
% %             for i = 1:10
% %                 plot([0.5],lnr0(:,1,1,i),'.b');
% %                 plot([2.0],lnr0(:,2,1,i),'.b');
% %                 plot([0.5],lnr0(:,1,2,i),'.r');
% %                 plot([2.0],lnr0(:,2,2,i),'.r');
% %             end
%             ylabel('ln(r_{45})','fontsize',18); %ylim([0 0.4]);
%             xlabel('Speed (rev/s)','fontsize',18);
%             legend('CW','CCW','Location','North'); grid on;
%             set(gca,'FontSize',20);
%             hold off;
            
%             outputt_ln0.P
%             outputt_ln45.P
            
        end
        
        function [] = get_lnrCoordinate_plots(this,subjectPlot,anovaPlot)
             
            speedVec = [1,2,3,1,2,3];
            dirVec = [1,1,1,2,2,2];
            
            speedSaveDex = [1,2,3,1,2,3];
            
            stiff = this.stiffDex;
            for subj = this.subjDex
                for speed = [1,2,3,4,5,6]
                    for trial = this.trialDex
                        
                        % Raw Dependent Measure
                        %this.test{subj,speed,trial}.majorAngle;
                        % Masure from bottom CCW direction
                        lnr0(trial,speedSaveDex(speed),dirVec(speed),subj) = this.test{subj,speed,trial}.lnrCoordinates(1);
                        lnr45(trial,speedSaveDex(speed),dirVec(speed),subj) = this.test{subj,speed,trial}.lnrCoordinates(2);
                        
                    end
                end
            end
            
            % Take out major outliers
            [lnr0] = this.remove_outliers(lnr0);
            [lnr45] = this.remove_outliers(lnr45);
            
            [PC1,this.outputt_lnr_PC] = this.get_lnEig(lnr0,lnr45,anovaPlot);
            this.get_lnr0lnr45(lnr0,lnr45,PC1,subjectPlot);
   
        end
        
        function [] = get_lnrCoordinates_plots_allStiff(this)
           
            outputt_lnr_PC = this.outputt_lnr_PC;
            
            % Plot ln0
            e1 = errorbar([0.075,0.5,2.0], outputt_lnr_PC.ave_CW, outputt_lnr_PC.dev_CW,'-','color',this.colorr{1},'linewidth',2.5,'markersize',30);
            set([e1.Bar,e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*(this.stiffDex/9)]);

            scatter([0.075,0.5,2.0], outputt_lnr_PC.ave_CW,50,'MarkerFaceColor',this.colorr{1},'MarkerEdgeColor',this.colorr{1},'MarkerEdgeAlpha',(this.stiffDex/9),'MarkerFaceAlpha',(this.stiffDex/9));


            
            e1 = errorbar([0.075,0.5,2.0], outputt_lnr_PC.ave_CCW, outputt_lnr_PC.dev_CCW,'-','color',this.colorr{2},'linewidth',2.5,'markersize',30);
            set([e1.Bar,e1.Line], 'ColorType', 'truecoloralpha', 'ColorData', [e1.Line.ColorData(1:3); 255*(this.stiffDex/9)]);

            scatter([0.075,0.5,2.0], outputt_lnr_PC.ave_CCW,50,'MarkerFaceColor',this.colorr{2},'MarkerEdgeColor',this.colorr{2},'MarkerEdgeAlpha',(this.stiffDex/9),'MarkerFaceAlpha',(this.stiffDex/9));

%             plot([0.075,0.5,2.0], outputt_lnr_PC.ave_CW,'.','color',this.colorr{1},'markersize',30);
%             plot([0.075,0.5,2.0], outputt_lnr_PC.ave_CCW,'.','color',this.colorr{2},'markersize',30);
                        ylim([-1.5 1.5]); %yticks(gca,[135:45:225]);
                        xlim([0 2.075]); xticks(gca,[0,0.5,2]);
            
            ylabel('ln(r)_{PC}','fontsize',18); %ylim([0 0.4]);
            xlabel('Speed (rev/s)','fontsize',18);
%             legend('CW','CCW','Location','North'); 
            set(gca,'FontSize',20,'linewidth',2.5);
            hold off;
            
        end
        
        function [] = get_CVplots(this)
            
            speedVec = [1,2,3,1,2,3];
            dirVec = [1,1,1,2,2,2];
            
            speedSaveDex = [1,2,3,1,2,3];
            
            
            skiDex = [7, 2, 3;...
                7, 2, 6;...
                7, 5, 14;...
                7, 5, 21;...
                9, 2, 7;...
                9, 2, 8;...
                9, 2, 11;...
                9, 5, 12;...
                10, 5, 14;...
                10, 6, 6];

            count = 1;
            stiff = this.stiffDex;
            for subj = this.subjDex
                for speed = [1,2,3,4,5,6]
%                     figure; title(speed); hold on;
                    clear v_bin v0_bin
                    for trial = this.trialDex
                        
                        if(sum(sum(skiDex == [subj, speed,trial],2)==3)~=1)
                        % Raw Dependent Measure
                        %this.test{subj,speed,trial}.majorAngle;
                        % Masure from bottom CCW direction
                        % Convert all to (m/s)
                        
                        lc = this.test{subj,speed,trial}.lc;
                        sfrq = this.test{subj,speed,trial}.sfrq;
                        d1 = this.test{subj,speed,trial}.d1;
                        d2 = this.test{subj,speed,trial}.d2;
                        
                        X = this.test{subj,speed,trial}.X;
                        X_0 = this.test{subj,speed,trial}.X_0;
                        
                        % Cut data to one cycle
                        thcp = wrapTo2Pi(atan2(X(:,2)-d2,X(:,1)-d1));
                        
%                         dexReset = find(abs(diff(thcp)) > 4);
%                          if(length(dexReset)<2)
%                              dexRange = [1:length(thcp)];
%                              problem(count,:) =  [subj,speed,trial];
%                              count = count + 1;
%                          else
%                              dexRange = [dexReset(end-1):dexReset(end)];
%                          end
%                          
%                          X = X(dexRange,:);
%                          X_0 = X_0(dexRange,:);
                        
                        % Cut initalization error
                        if(speedVec(speed) == 1)
                            thcp = thcp(75:end);
                            X = X(75:end,:);
                            X_0 = X_0(75:end,:);
                        end
                                                
%                         cf = 25; % cutoff freqnency
%                         [b,a] = butter(4,cf/(sfrq/2)); % make filter
%                         X(:,1) = filtfilt(b,a,X(:,1)); % apply fitler
%                         X(:,2) = filtfilt(b,a,X(:,2)); % apply fitler
%                         X_0(:,1) = filtfilt(b,a,X_0(:,1)); % apply fitler
%                         X_0(:,2) = filtfilt(b,a,X_0(:,2)); % apply fitler
            
                        X_dot = sfrq*diff(X);
                        X_dot(end+1,:) = X_dot(end,:);
                        X_dot_0 = sfrq*diff(X_0);
                        X_dot_0(end+1,:) = X_dot_0(end,:);

                        
                        v = sqrt((X_dot(:,1).^2 + X_dot(:,2).^2));
                        v0 = sqrt((X_dot_0(:,1).^2 + X_dot_0(:,2).^2));  
                        
                        [ X_bin, tmp, tmp_v_bin ] = this.customBinning(X(:,1), X(:,2), v, thcp);
                        [ X_bin, tmp, tmp_v0_bin ] = this.customBinning(X_0(:,1), X_0(:,2), v0, thcp);
                        
                        if ~exist('v_bin','var') % Catch first iteration and define v_bin
                            for i = 1:length(tmp_v_bin)
                                v_bin{i} = tmp_v_bin{i};
                                v0_bin{i} = tmp_v0_bin{i};
                            end
                        else
                            for i = 1:length(tmp_v_bin)
                                L1 = length(tmp_v_bin{i});
                                v_bin{i}(end+1:end+L1) = tmp_v_bin{i};
                                L2 = length(tmp_v0_bin{i});
                                v0_bin{i}(end+1:end+L2) = tmp_v0_bin{i};
                                clear L1 L2
                            end
                        end
                        
%                           v = this.test{subj,speed,trial}.thcv;
%                           v0 = this.test{subj,speed,trial}.Vt_0;
                          
%                           cv_v(trial,speedSaveDex(speed),dirVec(speed),subj) = std(v)/abs(mean(v));
%                           cv_v0(trial,speedSaveDex(speed),dirVec(speed),subj) = std(v0)/abs(mean(v0));

%                           cv_v(trial,speedSaveDex(speed),dirVec(speed),subj) = std(v)/abs(mean(v));
%                           cv_v0(trial,speedSaveDex(speed),dirVec(speed),subj) = std(v0)/abs(mean(v0));
                          
%                           subplot(2,1,1); plot(v_bin);hold on;
%                           subplot(2,1,2); plot(v0_bin);hold on;

%                         figure;
%                         subplot(2,1,1); %plot(v); hold on;
%                         plot(lc*this.test{subj,speed,trial}.thcv);hold on; pause; close all;
%                         
%                         subplot(2,1,2); 
%                         plot(v0); hold on;
% %                         plot(this.test{subj,speed,trial}.Vt_0);hold on;
%                         title([subj,speed,trial,abs(mean(v0)),std(v0),std(v0)/abs(mean(v0))]); %pause; close all
                        
                    end
                    % compute mean
                    
%                     figure; 
%                     for i = 1:length(v_bin)
%                         subplot(2,1,1); plot(i.*ones(size(v_bin{i})),v_bin{i},'.'); hold on;
%                         subplot(2,1,2); plot(i.*ones(size(v0_bin{i})),v0_bin{i},'.');hold on;
%                     end
                    
                    for i = 1:length(v_bin)
                        v_mean(i) = mean(v_bin{i});
                        v_std(i) = std(v_bin{i});
                        v0_mean(i) = mean(v0_bin{i});
                        v0_std(i) = std(v0_bin{i});
                    end
                    
                    cv_v(speedSaveDex(speed),dirVec(speed),subj) = mean(v_std./v_mean);
                    cv_v0(speedSaveDex(speed),dirVec(speed),subj) = mean(v0_std./v0_mean);
                    end
                end
            end
            
%             figure;
%             set(gca,'FontSize',20,'linewidth',2.5);hold on;
%             errorbar([0.075,0.5,2.0],mean(cv_v(:,1,:),3),std(cv_v(:,1,:),0,3),this.colorr(1),'linewidth',2.5);hold on;
%             errorbar([0.075,0.5,2.0],mean(cv_v(:,2,:),3),std(cv_v(:,2,:),0,3),this.colorr(2),'linewidth',2.5);
%             ylabel('CV of V','fontsize',18); ylim([0 0.4]);
%             xlabel('Speed (rev/s)','fontsize',18);
%             legend('CW','CCW','Location','North'); 
%             
%             figure;
%             set(gca,'FontSize',20,'linewidth',2.5);hold on;
%             errorbar([0.075,0.5,2.0],mean(cv_v0(:,1,:),3),std(cv_v0(:,1,:),0,3),this.colorr(1),'linewidth',2.5); hold on;
%             errorbar([0.075,0.5,2.0],mean(cv_v0(:,2,:),3),std(cv_v0(:,2,:),0,3),this.colorr(2),'linewidth',2.5);
%             ylabel('CV of V_0','fontsize',18); ylim([0 0.4]);
%             xlabel('Speed (rev/s)','fontsize',18);
%             legend('CW','CCW','Location','North'); 
            
            % Take out major outliers
%             [cv_v] = this.remove_outliers(cv_v);
%             [cv_v0] = this.remove_outliers(cv_v0);
            % Remove outlier function can not be used with binning approach
            % outliers have to be removed before binning by postion
            
            cv_v_tmp = cv_v;
            cv_v = [];
            cv_v(1,:,:,:) = cv_v_tmp;
            clear cv_v_tmp
            
            cv_v0_tmp = cv_v0;
            cv_v0 = [];
            cv_v0(1,:,:,:) = cv_v0_tmp;
            clear cv_v0_tmp
            
            %% Run ANOVA on PC coordinates
            [outputt_cv_v] = this.threeWayANOVAJH(cv_v);
            [outputt_cv_v0] = this.threeWayANOVAJH(cv_v0);
            
            
            %% ttest
%             size(cv_v)
            
            % Slow to medium
            [H,P,CI] = ttest(squeeze(cv_v(1,1,1,:)), squeeze(cv_v(1,2,1,:)));
            [H,P,CI] = ttest(squeeze(cv_v(1,1,2,:)), squeeze(cv_v(1,2,2,:)));

            % Medium to fast
            [H,P,CI] = ttest(squeeze(cv_v(1,2,1,:)), squeeze(cv_v(1,3,1,:)));
            [H,P,CI] = ttest(squeeze(cv_v(1,2,2,:)), squeeze(cv_v(1,3,2,:)));


            %% Plot outpue
            % Plot cv_v
            figure; hold on;
            errorbar([0.075,0.5,2.0], outputt_cv_v.ave_CW, outputt_cv_v.dev_CW,'.-','color',this.colorr{1},'linewidth',2.5,'markersize',30);
            errorbar([0.075,0.5,2.0], outputt_cv_v.ave_CCW,outputt_cv_v.dev_CCW,'.-','color',this.colorr{2},'linewidth',2.5,'markersize',30);
%             plot([0.075,0.5,2.0], outputt_cv_v.ave_CW,'.-b','markersize',30);
%             plot([0.075,0.5,2.0], outputt_cv_v.ave_CCW,'.-r','markersize',30);
                        %ylim([-0.7 0.7]); %yticks(gca,[135:45:225]);
                        %xlim([0 2.075]); xticks(gca,[0,0.5,2]);
            
%             for i = 1:10
%                 plot([0.075],mean(cv_v(:,1,1,i)),'ob');
%                 plot([0.5],  mean(cv_v(:,2,1,i)),'ob');
%                 plot([2.0],  mean(cv_v(:,3,1,i)),'ob');
%                 plot([0.075],mean(cv_v(:,1,2,i)),'ob');
%                 plot([0.5],  mean(cv_v(:,2,2,i)),'or');
%                 plot([2.0],  mean(cv_v(:,3,2,i)),'or');
%             end
            ylabel('CV of V','fontsize',18); ylim([0 0.4]); xticks([0,0.5,2]);
            xlabel('Speed (rev/s)','fontsize',18);
            legend('CW','CCW','Location','North'); 
            set(gca,'FontSize',20,'linewidth',2.5);
            hold off;
            
            % Plot cv_v0
            figure; hold on;
            errorbar([0.075,0.5,2.0], outputt_cv_v0.ave_CW, outputt_cv_v0.dev_CW,'color',this.colorr{1},'linewidth',2.5);
            errorbar([0.075,0.5,2.0], outputt_cv_v0.ave_CCW,outputt_cv_v0.dev_CCW,'color',this.colorr{2},'linewidth',2.5);
            plot([0.075,0.5,2.0], outputt_cv_v0.ave_CW,'.','color',this.colorr{1},'markersize',30);
            plot([0.075,0.5,2.0], outputt_cv_v0.ave_CCW,'.','color',this.colorr{2},'markersize',30);
%                         ylim([0 1]); %yticks(gca,[135:45:225]);
                        xlim([0 2.075]); xticks(gca,[0,0.5,2]);
            
            for i = 1:10
                plot([0.075],mean(cv_v0(:,1,1,i)),'o','color',this.colorr{1});
                plot([0.5],  mean(cv_v0(:,2,1,i)),'o','color',this.colorr{1});
                plot([2.0],  mean(cv_v0(:,3,1,i)),'o','color',this.colorr{1});
                plot([0.075],mean(cv_v0(:,1,2,i)),'o','color',this.colorr{2});
                plot([0.5],  mean(cv_v0(:,2,2,i)),'o','color',this.colorr{2});
                plot([2.0],  mean(cv_v0(:,3,2,i)),'o','color',this.colorr{2});
            end
            ylabel('CV of V_0','fontsize',18); ylim([0 0.4]);
            xlabel('Speed (rev/s)','fontsize',18);
            legend('CW','CCW','Location','North'); 
            set(gca,'FontSize',20,'linewidth',2.5);
            hold off;
            
            this.outputt_cv_v = outputt_cv_v;
                        
        end
                
        function [] = get_subjectAveZFTplots(this)
            % Filename:	get_subjectAveZFTplots.m
            % Author:  James Hermus
            % Date:		July 2 2019
            % Description:  This script makes the averaging figures for the curvature
            % velocity paper
            
            test = this.test;
            
            for subj = this.subjDex
                for speed = this.speedDex
                    X_0_tot{subj,speed} = [];
                    Vt_0_tot{subj,speed} = [];
                    thcp_0_tot{subj,speed} =[];
                    for trial = this.trialDex
                        X_0_tot{subj,speed} = [X_0_tot{subj,speed}; test{subj,speed,trial}.X_0];
                        Vt_0_tot{subj,speed} = [Vt_0_tot{subj,speed}; test{subj,speed,trial}.Vt_0];
                        thcp_0_tot{subj,speed} = [thcp_0_tot{subj,speed}; test{subj,speed,trial}.thcp_0];
                    end
                    
                    [ X_0_tot_bin{subj,speed}, Vt_0_tot_bin{subj,speed} ] = this.customBinning(X_0_tot{subj,speed}(300:end,1), X_0_tot{subj,speed}(300:end,2), Vt_0_tot{subj,speed}(300:end), thcp_0_tot{subj,speed}(300:end));
                    r_0_bin{subj,speed} = sqrt((X_0_tot_bin{subj,speed}(1,:) - test{subj,speed,trial}.d1).^2 + (X_0_tot_bin{subj,speed}(2,:) - test{subj,speed,trial}.d2).^2);
                    
                end
            end
            
            % Make plots
            LW = 2; % Line Width
            MS = 10; % Marker Size
            FS = 18; % Font Size
            
            theta = linspace(0,2*pi,length(r_0_bin{this.subjDex(1),this.speedDex(1)}));
            
            % Plot all subjects on same plot
            for speed = this.speedDex
                figure; hold on;
                for subj = this.subjDex
                    
                    pointsize = 50;
                    lc = 0.1029;
                    %         scatter( subjData(subj,speed).r_bin.*cos(theta),subjData(subj,speed).r_bin.*sin(theta), pointsize,subjData(subj,speed).vel_bin,'linewidth',LW); hcb = colorbar;
                    vel_max_min_scaled = (Vt_0_tot_bin{subj,speed} - min(Vt_0_tot_bin{subj,speed}))/(max(Vt_0_tot_bin{subj,speed})- min(Vt_0_tot_bin{subj,speed}));
                    scatter( r_0_bin{subj,speed}.*cos(theta),r_0_bin{subj,speed}.*sin(theta), pointsize,vel_max_min_scaled,'filled','linewidth',LW); hcb = colorbar;
                    plot(lc*cos(0:0.01:2*pi),lc*sin(0:0.01:2*pi),'--k','linewidth',2.5);
                    ylabel(hcb,'Normalized Speed');
                    ylabel('Y-Position (m)');
                    xlabel('X-Position (m)');axis equal; % grid on;
                    ylim([-0.18 0.18]);
                    xlim([-0.18 0.18]);
                    set(gca,'FontSize',28);
                    caxis([0 1]);
                end
                %         saveas(gcf, [folder,'TotNormZFT_',int2str(speed),'.eps']);
            end
            
        end
        
        function [outputt] = threeWayANOVAJH(this, Y )
            %Compute 3 Way Mixed ANOVA model:
            % A: speed (Fixed)
            % B: Direction (Fixed)
            % C: subject (Random)
            % Y_{i,j,k,l} = A + B + C + (AxB) + (AxC) + (BxC) + (AxBxC)
            
            [n,a,b,c] = size(Y);
            % n number of trials index i
            % a number of speeds index j
            % b number of directions index k
            % c number of subjects index l
            % Compute degrees of freedome (df)
            
            df.A = (a-1);
            df.B = (b-1);
            df.C = (c-1);
            df.A_B = (a-1)*(b-1);
            df.A_C = (a-1)*(c-1);
            df.B_C = (b-1)*(c-1);
            df.A_B_C = (a-1)*(b-1)*(c-1);
            df.Error = a*b*c*(n-1);
            df.T = n*a*b*c;
            
            % % Samity check Sum of all equals total - 1
            % df.A + df.B + df.C + df.A_B + df.A_C + df.B_C + df.A_B_C + df.Error
            % df.T - 1
            Y_bar = mean(mean(mean(mean(Y))));
            
            % Compute sums of squares (MS)
            SS.A = b*c*n*sum( ( mean(mean(mean(Y,1),3),4) - Y_bar ).^2 );
            SS.B = a*c*n*sum( ( mean(mean(mean(Y,1),2),4) - Y_bar ).^2 );
            SS.C = a*b*n*sum( ( mean(mean(mean(Y,1),2),3) - Y_bar ).^2 );
            SS.A_B = c*n*sum( sum( ( mean(mean(Y,1),4) - mean(mean(mean(Y,1),3),4) - mean(mean(mean(Y,1),2),4) + Y_bar ).^2 , 2 ), 3 );
            SS.A_C = b*n*sum( sum( ( mean(mean(Y,1),3) - mean(mean(mean(Y,1),3),4) - mean(mean(mean(Y,1),2),3) + Y_bar ).^2, 2 ), 4 );
            SS.B_C = a*n*sum( sum( ( mean(mean(Y,1),2) - mean(mean(mean(Y,1),2),4) - mean(mean(mean(Y,1),2),3) + Y_bar ).^2, 3 ), 4 );
            SS.A_B_C = n*sum(sum(sum( ( mean(mean(mean(Y,1),3),4) + mean(mean(mean(Y,1),2),4) + mean(mean(mean(Y,1),2),3) - ...
                mean(mean(Y,1),4) - mean(mean(Y,1),3) -mean(mean(Y,1),2) + mean(Y,1) - Y_bar ).^2 ,2),3),4);
            SS.Error = sum(sum(sum(sum( (Y - mean(Y,1) ).^2 ,1),2),3),4);
            SS.T = sum(sum(sum(sum( (Y - Y_bar ).^2 ,1),2),3),4);
            
            % % Samity check Sum of all equals total
            % SS.A + SS.B + SS.C + SS.A_B + SS.A_C + SS.B_C + SS.A_B_C + SS.Error
            % SS.T
            
            % Compute mean sums of squares
            MS.A = SS.A/df.A;
            MS.B = SS.B/df.B;
            MS.C = SS.C/df.C;
            MS.A_B = SS.A_B/df.A_B;
            MS.A_C = SS.A_C/df.A_C;
            MS.B_C = SS.B_C/df.B_C;
            MS.A_B_C = SS.A_B_C/df.A_B_C;
            MS.Error = SS.Error/df.Error;
            
            % Compute F ratio
            F.A = MS.A/MS.A_C;
            F.B = MS.B/MS.B_C;
            F.C = MS.C/MS.Error;
            F.A_B = MS.A_B/MS.A_B_C;
            F.A_C = MS.A_C/MS.Error;
            F.B_C = MS.B_C/MS.Error;
            F.A_B_C = MS.A_B_C/MS.Error;
            
            % Compute P value
            P.A = fpdf(F.A, df.A, df.A_C);
            P.B = fpdf(F.B, df.B, df.B_C);
            P.C = fpdf(F.C, df.C, df.Error);
            P.A_B = fpdf(F.A_B, df.A_B, df.A_B_C);
            P.A_C = fpdf(F.A_C, df.A_C, df.Error);
            P.B_C = fpdf(F.B_C, df.B_C, df.Error);
            P.A_B_C = fpdf(F.A_B_C, df.A_B_C, df.Error);
            
            outputt.df = df;
            outputt.SS = SS;
            outputt.MS = MS;
            outputt.F = F;
            outputt.P = P;
            
            % Plot info
            outputt.ave_CW = mean(mean(Y(:,:,1,:),1),4);
            outputt.dev_CW = std(squeeze(mean(Y(:,:,1,:),1))');
            outputt.ave_CCW = mean(mean(Y(:,:,2,:),1),4);
            outputt.dev_CCW = std(squeeze(mean(Y(:,:,2,:),1))');
            
        end
        
        % Bins the postion and exports the radius too for the ZFT plot
        function [ X_bin, vel_bin, vel_bin_tot ] = customBinning(this, x, y, vel, p)
            
            N = 200; % number of bins
            posEdges = linspace(0,2*pi,N+1);
            
            for j = 1:N
                
                binDex = find(p > posEdges(j) & p <= posEdges(j+1)); % Find index pos prosition with in each bin
                inda(j,1:length(binDex)) = binDex;                   % save index ast inda
                numPerBin(j) = length(binDex);
                vel_bin(j) = mean(vel(binDex));                  % take mean of velocity positon points with in each bin
                vel_bin_tot{j} = vel(binDex);
                
                % Find mean radius too
                x_bin(j) = mean(x(binDex));
                y_bin(j) = mean(y(binDex));
                r(j) = sqrt(mean(x(binDex)).^2+mean(y(binDex)).^2);
                
            end
            X_bin = [x_bin; y_bin];
            
        end
        
        function [] = multiColorLine(this,x,y,c,cmap)
            
            numPoints = numel(x);
            
            if nargin < 4
                cmap = parula;
            end
            cn = (c-min(c))/(max(c) - min(c));
            cn = ceil(cn*size(cmap,1));
            cn = max(cn,1);
            
            for i = 1:numPoints - 1
                line( x(i:i+1), y(i:i+1), 'color', cmap(cn(i),:), 'linewidth',5);
            end
            
        end
        
        function [depMeas] = remove_outliers(this,depMeas)
        
                    % Take out major outliers
            % Subject 5 outliers all over the place CCW medium
            % Subject 8 outliers all over the place CCW medium
            
            % Subject 7 outlier CW medium: 3, 6, CCW medium: 14, 21
            dexNoOutlier = find((this.trialDex ~= 3 )&( this.trialDex ~= 6));
            depMeas(3,1,1,7) = mean(squeeze(depMeas([dexNoOutlier],1,1,7)));
            depMeas(6,1,1,7) = mean(squeeze(depMeas([dexNoOutlier],1,1,7)));
            
            dexNoOutlier = find((this.trialDex ~= 14 )&( this.trialDex ~= 21));
            depMeas(14,1,2,7) = mean(squeeze(depMeas([dexNoOutlier],1,2,7)));
            depMeas(21,1,2,7) = mean(squeeze(depMeas([dexNoOutlier],1,2,7)));
            
            % Subject 9 outliers CW 7,8,11, CCW Medium 12
            dexNoOutlier = find((this.trialDex ~= 7 )&( this.trialDex ~= 8)&( this.trialDex ~= 11));
            depMeas(7,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            depMeas(8,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            depMeas(11,1,1,9) = mean(squeeze(depMeas([dexNoOutlier],1,1,9)));
            
            dexNoOutlier = find((this.trialDex ~= 12 ));
            depMeas(12,1,2,9) = mean(squeeze(depMeas([dexNoOutlier],1,2,9)));
            
            % Subject 10 outliers CCW medium 14, CCW fast 6
            dexNoOutlier = find((this.trialDex ~= 14 ));
            depMeas(14,1,2,10) = mean(squeeze(depMeas([dexNoOutlier],1,2,10)));
            
            dexNoOutlier = find((this.trialDex ~= 6 ));
            depMeas(6,2,2,10) = mean(squeeze(depMeas([dexNoOutlier],2,2,10)));
        
        end
        
    end
end

