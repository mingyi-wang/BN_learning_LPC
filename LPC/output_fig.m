figure 

subplot(2,3,1)
xlim([0 800])
n_measurements1=n_measurements;
n_measurements=100:100:700;
% s1=s;
s=100:100:700;
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
%set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
	hold all
	plot(n_measurements, AUC_PvsR_ggm_mean, '-d','color',[0 1 1])
% 	plot(n_measurements, AUC_PvsR_lowpc_mean, '-v', 'Color', [0 0.5 0])
    plot(n_measurements, AUC_PvsR_lowpc_oldcc_mean, '-+', 'Color', 'red')
    plot(n_measurements, AUC_PvsR_bn_mean, '-xb')
     plot(n_measurements, AUC_PvsR_rn_mean, '-o','Color',[1 0 1])
     
     plot(n_measurements, AUC_PvsR_pc_mean, '-v','color',[0 0.5 0.3])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_PvsR_o2pcorr_mean, '-^r')
% 	end
% 	plot(n_measurements, AUC_PvsR_ggm_mean, '-p', 'Color', [0 0.75 0.75])
% 	plot(n_measurements, AUC_PvsR_MI_mean, '-s', 'Color', [0.75 0 0.75])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, AUC_PvsR_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, AUC_PvsR_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
	xlabel('N. of measurements')
	ylabel('AUC(PvsR)')
	%ylim([0 0.7])
	title('(a)')
% 	if ~skip_lengthy_algorithms
% 		legend('R', 'R_{C1}', 'R_{C2}', 'R_{Call}', 'I', 'I_C', 'I_{DPI}', 'Location', 'NorthWest')
% 	else
% 		legend('R', 'R_{C1}', 'R_{Call}', 'I', 'I_{DPI}', 'Location', 'NorthWest')
% 	end
% 	if save_all
% 		saveas(gcf, [output_directory filenames_start '_AUC_PvsR.fig'])
% 	end
% 
 	
    subplot(2,3,4)
    xlim([0 800])
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
 	hold all
 	plot(n_measurements, TP_ggm_mean, '-d','color',[0 1 1])
% 	plot(n_measurements, TP_lowpc_mean, '-v', 'Color', [0 0.5 0])
    plot(n_measurements, TP_lowpc_oldcc_mean, '-+', 'Color', 'red')
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, TP_o2pcorr_mean, '-^r')
% 	end
 	plot(n_measurements, TP_bn_mean, '-xb')
    plot(n_measurements, TP_rn_mean, '-o','Color',[1 0 1])
 	plot(n_measurements, TP_pc_mean, '-v','color',[0 0.5 0.3])
% 	if ~skip_lengthy_algorithms
% 		plot(n_measurements, TP_minMI_mean, '-o', 'Color', [0.75 0.75 0])
% 	end
% 	plot(n_measurements, TP_DPI_mean, '-*', 'Color', [0.25 0.25 0.25])
%  	xlabel('N. of measurements')
%  	ylabel('TP for FP=20')
%  	title('(b)')
%     
     subplot(2,3,2);
     xlim([0 800])
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
 	hold all
 	errorbar(s, mean_num_shd, std_num_shd, '-+r')
 	errorbar(s, mean_num_shd1,std_num_shd1, '-v','color',[0 0.5 0.3])
   % errorbar(s, mean_num_shd2,std_num_shd2, '-xb')
 	xlabel('N. of measurements')
 	ylabel('Total number of structual errors')
 	title('(c)')
    
    subplot(2,3,3);
   xlim([0 800])
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
 	hold all
 	errorbar(s, mean_num_ea, std_num_ea,'-+r')
 	errorbar(s, mean_num_ea1,std_num_ea1, '-v','color',[0 0.5 0.3])
    %errorbar(s, mean_num_ea2,std_num_ea2,'-xb')
 	xlabel('N. of measurements')
 	ylabel('Number of extra arcs')
 	title('(d)')
    
    subplot(2,3,5);
    xlim([0 800])
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
 	hold all
 	errorbar(s, mean_num_ma, std_num_ma,'-+r')
 	errorbar(s, mean_num_ma1,std_num_ma1, '-v','color',[0 0.5 0.3])
    %errorbar(s, mean_num_ma2,std_num_ma2, '-xb')
 	xlabel('N. of measurements')
 	ylabel('Number of missing arcs')
 	title('(e)')
    
    subplot(2,3,6);
    xlim([0 800])
set(gca,'XTick',0:100:800 ,'XTickLabel',{'','10','20','50','100','200','500','1000',''},'fontsize',12); 
 	hold all
 	errorbar(s, mean_num_wo,std_num_wo, '-+r')
 	errorbar(s, mean_num_wo1, std_num_wo1,'-v','color',[0 0.5 0.3])
    %errorbar(s, mean_num_mo2, std_num_mo2,'-xb')
 	xlabel('N. of measurements')
 	ylabel('Number of arcs with reversed or no directions')
 	title(['(f)'])