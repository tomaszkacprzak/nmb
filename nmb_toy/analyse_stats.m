clean
runid='001002'
cd(runid)
pwd 

fontsize_legend = 10;
markersize_true_si = 20;
do_noise_repetition = 0;

% make directory for figures if not existent
!mkdir figures

% matlab iterates from 1, so add 1 to all python column numbers
cols=addToFields(ReadYaml('tables_info.yaml'),1);

% figure file format
figure_file_format = 'png';

% true sersic indices
true_si_range = [1,2,4]
for true_si = true_si_range

	% get the stats from different noise maps
	S=load('nmb_toy.stats.d.txt');

	true_e = 0.2;
	select = S(:,cols.stats.true_si)==true_si;
	S=S(select,:);

	cols.stats.model_si
	si = S(:,cols.stats.model_si);
	model_g1 = S(:,cols.stats.model_g1);
	noise_bfit_g1     = S(:,cols.stats.noise_bfit_g1); 
	noise_bfit_g1_stdm = S(:,cols.stats.noise_bfit_g1_stdm); 

	noise_real_g1     = S(:,cols.stats.noise_real_g1); 
	noise_real_g1_stdm = S(:,cols.stats.noise_real_g1_stdm); 

	% get the diff from same noise maps
	S=load('nmb_toy.stats.s.txt');
	select = S(:,cols.stats.true_si)==true_si;
	S=S(select,:);

	diff_g1 = S(:,cols.stats.diff_g1); 
	diff_g1_stdm = S(:,cols.stats.diff_g1_stdm); 
	diff_g1_stdv = S(:,cols.stats.diff_g1_stdv); 

	% FIGURE 1 --------------------------------------------------------------------------
	figure(1); clf
	plot(si,model_g1,'r.-');;hold on
	% plot the star for the true point
	plot(true_si,model_g1(si==true_si),'pc','MarkerSize',markersize_true_si)

	hold on
	errorbar(si,noise_real_g1,noise_real_g1_stdm,noise_real_g1_stdm,'b+-')
	errorbar(si,noise_bfit_g1,noise_bfit_g1_stdm,noise_bfit_g1_stdm,'mx-')
	if do_noise_repetition == 1
	errorbar(si,noise_bfit_g1+diff_g1,diff_g1_stdm+noise_bfit_g1_stdm,diff_g1_stdm+noise_bfit_g1_stdm,'gx-')
	end
	% errorbar(si,noise_bfit_g1+diff_g1,diff_g1_stdm,diff_g1_stdm,'cx-')
	% errorbar(si,noise_bfit_g1+diff_g1,diff_g1_stdv,diff_g1_stdv,'cx-')
	plot(get(gca(),'xlim'),ones(2,1)*true_e,'k')
	grid on

	ylabel('e')
	xlabel('si_j')
	% ylim([0.2 0.32])
	% ylim([0.22 0.4])
	legend('true g1')
	legend1 = ['(1 model bias) ' char(10) ' shear of true galaxy with si_t=0.2 when fitted with model with si_j'];
	legend11= ['true Sersic index'];
	legend2 = ['(2 total bias) ' char(10) ' mean shear on true galaxy with si_t=0.2 fitted with si_j'];
	legend3 = ['(3 noise bias on the bestfit + model bias) ' char(10) ' mean shear when true galaxy is a best fit with si_j to noiseless galaxy with si_t=0.2, noisy gals fitted with si_j'];
	legend4 = ['(3 + [2-3 from noise repetition]) ' char(10) ' difference between [total bias] and [noise bias on the best fit + model bias] was calculated using same noise maps and 5000 noise realisations'];
	if do_noise_repetition == 1
		legends = {legend1,legend11,legend2,legend3,legend4};
	else
		legends = {legend1,legend11,legend2,legend3};
	end
	hl=legend(legends,'Location','SouthOutside');
	
	set(hl,'FontSize',fontsize_legend)

	title(sprintf('shear, true sersic-index=%2.2f, true-e1=%2.2f, points outside are missing results',true_si,true_e))

	paper_position = [0 0 20 15];
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', paper_position);
	filename_fig = sprintf('figures/fig1_panel1_SItrue%d.%s',true_si,figure_file_format);
	print(filename_fig,sprintf('-d%s',figure_file_format))
	fprintf('saved %s\n',filename_fig)

	% FIGURE 2 --------------------------------------------------------------------------
	figure(2); clf
	frac_diff = (model_g1 - true_e)/true_e;
	plot(si,frac_diff,'r.-');hold on
	plot(true_si,frac_diff(si==true_si),'pc','MarkerSize',30)


	frac_diff = (noise_real_g1 - true_e)/true_e;
	frac_diff_std = noise_real_g1_stdm;
	errorbar(si,frac_diff,frac_diff_std,frac_diff_std,'b.-');hold on

	frac_diff = (noise_bfit_g1 - true_e)/true_e;
	frac_diff_std = noise_bfit_g1_stdm;
	errorbar(si,frac_diff,frac_diff_std,frac_diff_std,'m.-');hold on

	frac_diff = (noise_bfit_g1+diff_g1-true_e)/true_e;
	frac_diff_var = (noise_bfit_g1_stdm.^2+diff_g1_stdm.^2)./true_e.^2;
	frac_diff_std = sqrt(frac_diff_var);
	if do_noise_repetition == 1
	errorbar(si,frac_diff,frac_diff_std,frac_diff_std,'g.-');hold on
	end
	grid on

	xlabel('si_j')
	ylabel('(fractional bias (g_a - g_b)/g_b')

	legend1 = ['(1 model bias) ' char(10) ' shear of true galaxy with si_t=0.2 when fitted with model with si_j'];
	legend11= ['true Sersic index'];
	legend2 = ['(2 total bias) ' char(10) ' mean shear on true galaxy with si_t=0.2 fitted with si_j'];
	legend3 = ['(3 noise bias on the bestfit + model bias) ' char(10) ' mean shear when true galaxy is a best fit with si_j to noiseless galaxy with si_t=0.2, noisy gals fitted with si_j'];
	legend4 = ['(3 + [2-3 from noise repetition]) ' char(10) ' difference between [total bias] and [noise bias on the best fit + model bias] was calculated using same noise maps and 5000 noise realisations'];
	if do_noise_repetition == 1
		legends = {legend1,legend11,legend2,legend3,legend4};
	else
		legends = {legend1,legend11,legend2,legend3};
	end
	hl=legend(legends,'Location','SouthOutside');
	set(hl,'FontSize',fontsize_legend)
	title('fractional bias (g_a-g_b)/g_a')
	% ylim([-0.03 0.03])
	plot(get(gca(),'xlim'),ones(2,1)*0,'k')

	paper_position = [0 0 20 15];
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', paper_position);
	filename_fig = sprintf('figures/fig1_panel2_SItrue%d.%s',true_si,figure_file_format);
	print(filename_fig,sprintf('-d%s',figure_file_format))
	fprintf('saved %s\n',filename_fig)

	% FIGURE 3 --------------------------------------------------------------------------

	figure(3); clf
	frac_diff_diffnoise = (noise_real_g1 -  noise_bfit_g1)./true_e;
	frac_diff_diffnoise_var = (noise_real_g1_stdm.^2 + noise_bfit_g1_stdm.^2)/true_e.^2;
	frac_diff_diffnoise_std = sqrt(frac_diff_var);

	frac_diff_samenoise = diff_g1./true_e;
	frac_diff_samenoise_stdm = diff_g1_stdm./true_e;
	frac_diff_samenoise_stdv = diff_g1_stdv./true_e;

	s1=errorbar(si,frac_diff_diffnoise,frac_diff_diffnoise_std,frac_diff_diffnoise_std,'b.-');hold on
	s2=errorbar(si,frac_diff_samenoise,frac_diff_samenoise_stdm,frac_diff_samenoise_stdm,'r.-');hold on
	% ylim([-0.03 0.03])
	% plot the star for the true point
	s3=plot(true_si,0,'pc','MarkerSize',markersize_true_si)
	s4=plot(get(gca(),'xlim'),ones(2,1)*0,'k')
	grid on
	
	legend1 = ['[   (2 total bias)  - (3 noise bias on the bestfit + model bias) ]/e_{true}'];
	legend2 = ['[   (2 total bias)  - (3 noise bias on the bestfit + model bias) ]/e_{true} from same noise maps'];
	legend3 = ['true Sersic index'];
	legends = {legend1,legend2,legend3};
	series = {s1,s2,s3,s4};
    	hl=legend([s1,s2,s3,s4],legends,'Location','SouthOutside');
	set(hl,'FontSize',fontsize_legend)
	title('contribution to the fractional shear bias from noise and model bias interaction')

	paper_position = [0 0 20 15];
	set(gcf, 'PaperUnits', 'centimeters');
	set(gcf, 'PaperPosition', paper_position);
	filename_fig = sprintf('figures/fig1_panel3_SItrue%d.%s',true_si,figure_file_format);
	print(filename_fig,sprintf('-d%s',figure_file_format))
	fprintf('saved %s\n',filename_fig)


end

cd('..')