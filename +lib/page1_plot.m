function f = page1_plot(dwmri_file,bval_file,mask_file,stats_mask_file,sigma_est,outliers_file,chi_squared_plot,dwi_model_file,movement_params_file)
    % This does not do any reorientation of nifti for chi squared analysis
    % because output should be WRT storage order. I did this because
    % generally storage order is the order acquired by the scanner, and the
    % chi^2 is done slice-wise. 
    
    % Load data 
    dwmri_hdr = load_untouch_header_only(dwmri_file.get_path());
    dwmri_vol = nifti_utils.load_untouch_nii4D_vol_scaled(dwmri_file.get_path(),'double');
    bvals = bval_file.dlmread();    
    mask_vol = nifti_utils.load_untouch_nii_vol_scaled(mask_file.get_path(),'logical');    
    stats_mask_vol = nifti_utils.load_untouch_nii_vol_scaled(stats_mask_file.get_path(),'logical');    
    outliers_vol = nifti_utils.load_untouch_nii_vol_scaled(outliers_file.get_path(),'logical'); 
    dwi_model_vol = nifti_utils.load_untouch_nii4D_vol_scaled(dwi_model_file.get_path(),'double');   
    xform_RAS = nifti_utils.get_voxel_RAS_xform(dwmri_file.get_path());
    
    % Get additional data
    dwi_vol = dwmri_vol(:,:,:,bvals ~= 0);    
    b0_vol = nanmean(dwmri_vol(:,:,:,bvals == 0),4);        
    
    % Get dwmri base plot ------------------------------------------------%
    split_file = strsplit(mfilename('fullpath'),filesep); % Used to get version, should be safe
    [f,pos_header,pos_info,pos_footer] = dwmri_info_plot({'DTI STATS', ...
                                                         ['   Scan Description: ' dwmri_hdr.hist.descrip], ...
                                                         ['   Gradient Directions: ' num2str(length(find(bvals)))], ...
                                                         ['   Slice Dimensions: ' num2str(dwmri_hdr.dime.dim(2)) ',' num2str(dwmri_hdr.dime.dim(3))], ...
                                                         ['   Slices: ' num2str(dwmri_hdr.dime.dim(4))], ...
                                                         ['   Voxel Resolution: ' num2str(round(dwmri_hdr.dime.pixdim(2),1)) ' x ' num2str(round(dwmri_hdr.dime.pixdim(3),1)) ' x ' num2str(round(dwmri_hdr.dime.pixdim(4),1))], ...
                                                         'Justin Blaber (based on code from Carolyn B. Lauzon', ...
                                                         'justin.blaber@vanderbilt.edu', ...
                                                         'Simultaneous Analysis and Quality Assurance for Diffusion Tensor Imaging (Lauzon et al 2013)', ...
                                                         'bennett.landman@vanderbilt.edu', ...
                                                         split_file{end-2}}); %#ok<ASGLU>
                   
    % Set up axes --------------------------------------------------------%
    padding = 0.01;
    chi_squared_area_height = 0.50*(pos_info(2)-(pos_footer(2)+pos_footer(4))-2*padding);
    chi_squared_area_width = 1-2*padding;
    chi_squared_header_height = 0.011;
    chi_squared_fits_height = chi_squared_area_height - 3*padding - 2*chi_squared_header_height;
    chi_squared_fits_width = 0.175*chi_squared_area_width;    
    chi_squared_plot_height = 0.840*(chi_squared_fits_height);
    chi_squared_plot_width = 0.475*(chi_squared_area_width-2*chi_squared_fits_width);    
    chi_squared_plot_left_width = 0.50*((chi_squared_area_width - 2*chi_squared_fits_width - chi_squared_plot_width)/2);
    chi_squared_plot_top_height = padding + 2*chi_squared_header_height + (chi_squared_fits_height - chi_squared_plot_height)/2;
    
    % Create chi^2 axes
    pos_chi_squared_fits_worst = [padding pos_footer(2)+pos_footer(4)+padding chi_squared_fits_width chi_squared_fits_height];
    axes_chi_squared_fits_worst = axes('Position',pos_chi_squared_fits_worst,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    pos_chi_squared_fits_best = [padding+chi_squared_area_width-chi_squared_fits_width pos_footer(2)+pos_footer(4)+padding chi_squared_fits_width chi_squared_fits_height];
    axes_chi_squared_fits_best = axes('Position',pos_chi_squared_fits_best,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    pos_chi_squared_plot = [padding+(chi_squared_area_width-chi_squared_plot_width)/2 pos_footer(2)+pos_footer(4)+padding+(chi_squared_fits_height-chi_squared_plot_height)/2 chi_squared_plot_width chi_squared_plot_height];
    axes_chi_squared_plot = axes('Position',pos_chi_squared_plot,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);    
    pos_chi_squared_plot_left = [pos_chi_squared_plot(1)-chi_squared_plot_left_width-padding pos_chi_squared_plot(2) chi_squared_plot_left_width pos_chi_squared_plot(4)];
    axes_chi_squared_plot_left = axes('Position',pos_chi_squared_plot_left,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    pos_chi_squared_plot_top = [pos_chi_squared_plot(1) pos_chi_squared_plot(2)+pos_chi_squared_plot(4)+padding pos_chi_squared_plot(3) chi_squared_plot_top_height];
    axes_chi_squared_plot_top = axes('Position',pos_chi_squared_plot_top,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);    
    
    % Create outliers axis
    pos_outliers_plot_top = [pos_chi_squared_plot_top(1) pos_footer(2)+pos_footer(4)+2*padding+chi_squared_area_height pos_chi_squared_plot_top(3) pos_chi_squared_plot_top(4)];
    axes_outliers_plot_top = axes('Position',pos_outliers_plot_top,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    
    % Make a plot of b-value vs median signal intensity. This is helpful
    % for debugging if b-values are mixed up wrt dwi - a common error.
    pos_bval_intensity_plot = [pos_outliers_plot_top(1) pos_outliers_plot_top(2)+pos_outliers_plot_top(4)+3*padding pos_outliers_plot_top(3) pos_outliers_plot_top(4)];
    axes_bval_intensity_plot_top = axes('Position',pos_bval_intensity_plot,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
       
    % Translation and rotation axes
    pos_plot_trans = [pos_bval_intensity_plot(1) pos_bval_intensity_plot(2)+pos_bval_intensity_plot(4)+3*padding pos_bval_intensity_plot(3) pos_bval_intensity_plot(4)];
    axes_plot_trans = axes('Position',pos_plot_trans,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    pos_plot_rot = [pos_plot_trans(1) pos_plot_trans(2)+pos_plot_trans(4)+3*padding pos_plot_trans(3) pos_plot_trans(4)];
    axes_plot_rot = axes('Position',pos_plot_rot,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    
    % Mask axes
    mask_header = 0.011;
    mask_area_height = pos_info(2)-(pos_footer(2)+pos_footer(4)+padding+chi_squared_area_height)-5*padding;      
    mask_height = (mask_area_height-2*padding)/3;
    mask_width = 0.40*(chi_squared_area_width-chi_squared_plot_width)/2;
    
    pos_mask_sagittal = [(padding+(chi_squared_area_width-chi_squared_plot_width)/2-mask_width)/2 pos_footer(2)+pos_footer(4)+2*padding+chi_squared_area_height mask_width mask_height];
    axes_mask_sagittal = axes('Position',pos_mask_sagittal,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    
    pos_mask_coronal = [pos_mask_sagittal(1) pos_mask_sagittal(2)+pos_mask_sagittal(4)+padding mask_width mask_height];
    axes_mask_coronal = axes('Position',pos_mask_coronal,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
    
    pos_mask_axial = [pos_mask_coronal(1) pos_mask_coronal(2)+pos_mask_coronal(4)+padding mask_width mask_height];
    axes_mask_axial = axes('Position',pos_mask_axial,'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);
       
    % These axes are used for plotting stuff between other axes. You have
    % to make an axes over everything to accomplish this.
    axes_overlay = axes('Position',[0 0 1 1],'Color','none','Xlim',[0 1],'Ylim',[0 1],'xtick',[],'ytick',[],'XColor','w','YColor','w','Parent',f);    
    hold(axes_overlay,'on'); % Set hold here
        
    % End of axes --------------------------------------------------------%
    
    % chi^2 headers
    
    % Worst fits
    worst_model_pos = [pos_chi_squared_fits_worst(1)+padding pos_chi_squared_fits_worst(2)+pos_chi_squared_fits_worst(4)+padding (pos_chi_squared_fits_worst(3)-3*padding)/2 chi_squared_header_height];
    uicontrol('style','text','units','normalized','String',{'Model'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',worst_model_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    worst_data_pos = [worst_model_pos(1)+worst_model_pos(3)+padding worst_model_pos(2) worst_model_pos(3) worst_model_pos(4)];
    uicontrol('style','text','units','normalized','String',{'Data'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',worst_data_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    worst_fits_pos = [worst_model_pos(1) worst_model_pos(2)+worst_model_pos(4)+padding 2*worst_model_pos(3)+padding worst_model_pos(4)];
    uicontrol('style','text','units','normalized','String',{'Worst Fits'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',worst_fits_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    % Best fits
    best_data_pos = [pos_chi_squared_fits_best(1)+padding pos_chi_squared_fits_best(2)+pos_chi_squared_fits_best(4)+padding (pos_chi_squared_fits_best(3)-3*padding)/2 chi_squared_header_height];
    uicontrol('style','text','units','normalized','String',{'Data'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',best_data_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    best_model_pos = [best_data_pos(1)+best_data_pos(3)+padding best_data_pos(2) best_data_pos(3) best_data_pos(4)];
    uicontrol('style','text','units','normalized','String',{'Model'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',best_model_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    best_fits_pos = [best_data_pos(1) best_data_pos(2)+best_data_pos(4)+padding 2*best_data_pos(3)+padding best_data_pos(4)];
    uicontrol('style','text','units','normalized','String',{'Best Fits'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',best_fits_pos, ...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
            
    % Do chi^2 plot
    imagesc(chi_squared_plot,'parent',axes_chi_squared_plot,[0 0.2]);
    set(axes_chi_squared_plot,'Ydir','Normal','YAxisLocation','right');
    colormap(axes_chi_squared_plot,'jet');    
    xlabel(axes_chi_squared_plot,'Gradient #','FontUnits','normalized','fontsize',1/25,'fontweight','bold');
    ylabel(axes_chi_squared_plot,'Slice #','FontUnits','normalized','fontsize',1/25,'fontweight','bold');
    set(axes_chi_squared_plot,'FontSize',8);
    % Colorbar - this will resize axes; so you must set the size again
    c = colorbar(axes_chi_squared_plot);
    set(axes_chi_squared_plot,'Position',pos_chi_squared_plot);
    % Reset position of colorbar
    colorbar_pos = get(c,'Position');
    colorbar_pos(1) = pos_chi_squared_plot(1)+pos_chi_squared_plot(3)+7*padding;
    colorbar_pos(3) = 0.5*colorbar_pos(3);
    set(c,'Position',colorbar_pos,'YTick',[0 0.2]);
    drawnow

    % Do median plot left of chi^2
    plot(axes_chi_squared_plot_left,nanmedian(chi_squared_plot,2),'LineWidth',1,'color','black');
    view(axes_chi_squared_plot_left,-90,90); % Flip x and y axes
    set(axes_chi_squared_plot_left,'Xlim',[1 size(chi_squared_plot,1)],'Ylim',[0 0.2],'box','off','XAxisLocation','top','xtick',[],'ytick',[0 0.2]);
    ylabel(axes_chi_squared_plot_left,'chi^{2} (slice)','FontUnits','normalized','fontsize',1/25,'fontweight','bold');
    set(axes_chi_squared_plot_left,'FontSize',8);
    drawnow    
    
    % Do median plot on top of chi^2
    plot(axes_chi_squared_plot_top,nanmedian(chi_squared_plot,1),'LineWidth',1,'color','black');
    set(axes_chi_squared_plot_top,'Xlim',[1 size(chi_squared_plot,2)],'Ylim',[0 0.2],'box','off','YAxisLocation','right','xtick',[],'ytick',[0 0.2]);
    ylabel(axes_chi_squared_plot_top,'chi^{2} (grad)','FontUnits','normalized','fontsize',1/6,'fontweight','bold');
    set(axes_chi_squared_plot_top,'FontSize',8);
    drawnow
    
    % Get best and worst fit montages
    % Form a montage
    num_fits = 5;
    rows_slice = size(dwi_vol,2); % Slices are rot90'd when plotted so #rows and #cols are swapped
    cols_slice = size(dwi_vol,1);
    montage_worst = zeros(num_fits*rows_slice,2*cols_slice);
    montage_best = zeros(num_fits*rows_slice,2*cols_slice);
    % Get the best fits within 5th and 95th percentile of axial slices.
    % This prevents finding fits with the first/last slices, which are 
    % usually bad.
    num_slices = size(dwi_vol,3);
    slice_upperbound = round(0.95*(num_slices-1)) + 1;
    slice_lowerbound = round(0.05*(num_slices-1)) + 1;
    slice_range_idxs = round(linspace(slice_lowerbound,slice_upperbound,num_fits+1));
    % Store results
    worst_slices_idx = zeros(num_fits,1);
    worst_slices_grad_idx = zeros(num_fits,1);
    best_slices_idx = zeros(num_fits,1);
    best_slices_grad_idx = zeros(num_fits,1);
    max_chi_squared = zeros(num_fits,1);
    min_chi_squared = zeros(num_fits,1);
    % data used for plotting maximum and minimum values - grab both data
    % and model values
    worst_data = [];
    best_data = [];
    for i = 1:num_fits
        % Get slice with the worst/best chi^2 between slice_idxs(i) and 
        % slice_idxs(i+1)
        chi_squared_sub_plot = chi_squared_plot(slice_range_idxs(i):slice_range_idxs(i+1)-1,:); % Subtract 1 from last index so sub_plots do not overlap
        max_chi_squared(i) = max(chi_squared_sub_plot(:));
        min_chi_squared(i) = min(chi_squared_sub_plot(:));
        [worst_slices_idx(i),worst_slices_grad_idx(i)] = find(chi_squared_sub_plot == max_chi_squared(i),1);
        worst_slices_idx(i) = worst_slices_idx(i) + slice_range_idxs(i)-1; % Adjust it
        [best_slices_idx(i),best_slices_grad_idx(i)] = find(chi_squared_sub_plot == min_chi_squared(i),1);
        best_slices_idx(i) = best_slices_idx(i) + slice_range_idxs(i)-1; % Adjust it
        
        % Get slices
        worst_slice = rot90(dwi_vol(:,:,worst_slices_idx(i),worst_slices_grad_idx(i)));
        worst_model_slice = rot90(dwi_model_vol(:,:,worst_slices_idx(i),worst_slices_grad_idx(i)));
        best_slice = rot90(dwi_vol(:,:,best_slices_idx(i),best_slices_grad_idx(i)));
        best_model_slice = rot90(dwi_model_vol(:,:,best_slices_idx(i),best_slices_grad_idx(i)));
        
        % Store data
        worst_data = vertcat(worst_data,worst_slice(rot90(mask_vol(:,:,worst_slices_idx(i)))),worst_model_slice(rot90(mask_vol(:,:,worst_slices_idx(i))))); %#ok<AGROW>
        best_data = vertcat(best_data,best_slice(rot90(mask_vol(:,:,best_slices_idx(i)))),best_model_slice(rot90(mask_vol(:,:,best_slices_idx(i))))); %#ok<AGROW>
        
        % Store in montage - do this in reverse order since y-axis is
        % flipped
        montage_worst((num_fits-i)*rows_slice+1:(num_fits-i+1)*rows_slice,1:cols_slice) = worst_model_slice;
        montage_worst((num_fits-i)*rows_slice+1:(num_fits-i+1)*rows_slice,cols_slice+1:2*cols_slice) = worst_slice;
        montage_best((num_fits-i)*rows_slice+1:(num_fits-i+1)*rows_slice,1:cols_slice) = best_slice;
        montage_best((num_fits-i)*rows_slice+1:(num_fits-i+1)*rows_slice,cols_slice+1:2*cols_slice) = best_model_slice;
    end    
    
    % Plot montages
    imagesc(montage_worst,'parent',axes_chi_squared_fits_worst,[0 max(3*nanmedian(worst_data),3*nanmedian(best_data))]);
    colormap(axes_chi_squared_fits_worst,'gray');
    axis(axes_chi_squared_fits_worst,'off','image');
    hold(axes_chi_squared_fits_worst,'on');
    
    imagesc(montage_best,'parent',axes_chi_squared_fits_best,[0 max(3*nanmedian(worst_data),3*nanmedian(best_data))]);
    colormap(axes_chi_squared_fits_best,'gray');
    axis(axes_chi_squared_fits_best,'off','image');
    hold(axes_chi_squared_fits_best,'on');
    
    % Plot text over montages
    for i = 1:num_fits         
        % Worst fit text
        text(2*cols_slice,(num_fits-i+1)*rows_slice,num2str(round(max_chi_squared(i),2)),'Color','w','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6,'Parent',axes_chi_squared_fits_worst);         
        text(cols_slice,(num_fits-i+1)*rows_slice,['(' num2str(worst_slices_idx(i)) ',' num2str(worst_slices_grad_idx(i)) ')'],'Color','w','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6,'Parent',axes_chi_squared_fits_worst);   
        
        % Best fit text          
        text(cols_slice,(num_fits-i+1)*rows_slice,num2str(round(min_chi_squared(i),2)),'Color','w','HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',6,'Parent',axes_chi_squared_fits_best);         
        text(1,(num_fits-i+1)*rows_slice,['(' num2str(best_slices_idx(i)) ',' num2str(best_slices_grad_idx(i)) ')'],'Color','w','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6,'Parent',axes_chi_squared_fits_best);   
    end
    
    % Plot white lines and text over montages
    % Horizontal lines
    for i = 1:num_fits-1     
        line([1 2*cols_slice],[i*rows_slice i*rows_slice],'color','w','Parent',axes_chi_squared_fits_worst);     
        line([1 2*cols_slice],[i*rows_slice i*rows_slice],'color','w','Parent',axes_chi_squared_fits_best);       
    end
    % Vertical lines
    line([cols_slice cols_slice],[1 num_fits*rows_slice],'color','w','Parent',axes_chi_squared_fits_worst);     
    line([cols_slice cols_slice],[1 num_fits*rows_slice],'color','w','Parent',axes_chi_squared_fits_best);    
    drawnow
        
    % Plot lines connecting values in chi^2 plot to montages; must update
    % position since axes are "image'd"
    pos_chi_squared_fits_worst = plotboxpos(axes_chi_squared_fits_worst);
    pos_chi_squared_fits_best = plotboxpos(axes_chi_squared_fits_best);
    line_offset = 0.0075;
    for i = 1:num_fits
        % Plot worst points
        worst_plot_x1 = (worst_slices_grad_idx(i)-0.5)/size(chi_squared_plot,2)*pos_chi_squared_plot(3)+pos_chi_squared_plot(1);
        worst_plot_y1 = (worst_slices_idx(i)-0.5)/size(chi_squared_plot,1)*pos_chi_squared_plot(4)+pos_chi_squared_plot(2);
                
        plot(axes_overlay,worst_plot_x1,worst_plot_y1,'ks','MarkerSize',12,'LineWidth',0.5);
        
        % Plot line between worst points and corresponding data slice
        worst_plot_x2 = pos_chi_squared_fits_worst(1) + pos_chi_squared_fits_worst(3);
        worst_plot_y2 = pos_chi_squared_fits_worst(2) + (i-0.5)/num_fits*pos_chi_squared_fits_worst(4);
        
        line([worst_plot_x1-line_offset worst_plot_x2+line_offset/2],[worst_plot_y1 worst_plot_y2],'color',[0 0 0 0.2],'Parent',axes_overlay);
        
        % Plot best points
        best_plot_x1 = (best_slices_grad_idx(i)-0.5)/size(chi_squared_plot,2)*pos_chi_squared_plot(3)+pos_chi_squared_plot(1);
        best_plot_y1 = (best_slices_idx(i)-0.5)/size(chi_squared_plot,1)*pos_chi_squared_plot(4)+pos_chi_squared_plot(2);
                
        plot(axes_overlay,best_plot_x1,best_plot_y1,'ks','MarkerSize',12,'LineWidth',0.5);  
        
        % Plot line between best points and corresponding data slice      
        best_plot_x2 = pos_chi_squared_fits_best(1);
        best_plot_y2 = pos_chi_squared_fits_best(2) + (i-0.5)/num_fits*pos_chi_squared_fits_best(4);
        
        line([best_plot_x1+line_offset best_plot_x2-line_offset/2],[best_plot_y1 best_plot_y2],'color',[0 0 0 0.2],'Parent',axes_overlay);
    end
    
    % Plot outliers above chi^2 gradient histogram
    outliers_grad_pct = squeeze(sum(sum(sum(outliers_vol,1),2),3))./sum(mask_vol(:)) * 100;    
    ylim_outliers_min = 0;
    ylim_outliers_max = ceil(max(outliers_grad_pct(:))*10)/10; % ceils to next decimal
    plot(axes_outliers_plot_top,outliers_grad_pct,'LineWidth',1,'color','black');
    set(axes_outliers_plot_top,'Xlim',[1 length(outliers_grad_pct)],'Ylim',[ylim_outliers_min ylim_outliers_max],'box','off','YAxisLocation','right','xtick',[],'ytick',[ylim_outliers_min ylim_outliers_max],'yticklabel',{[num2str(ylim_outliers_min) ' %'],[num2str(ylim_outliers_max) ' %']},'FontUnits','normalized','fontsize',1/6);       
    title(axes_outliers_plot_top,'Outliers');
    
    % Plot bvalues and median signal intensity
    median_intensities = zeros(size(dwi_vol,4),1);    
    for i = 1:size(dwi_vol,4)
        dwmri_vol = dwi_vol(:,:,:,i);
        median_intensities(i) = nanmedian(dwmri_vol(mask_vol));
    end
    
    [hAx,l1,l2] = plotyy(axes_bval_intensity_plot_top,1:size(dwi_vol,4),bvals(bvals ~= 0),1:size(dwi_vol,4),median_intensities);
    set(l1,'LineWidth',1);
    set(l2,'LineWidth',1);
    set(hAx(1),'Xlim',[1 size(dwi_vol,4)],'box','off','xtick',[],'ytick',get(hAx(1),'YLim'),'FontUnits','normalized','fontsize',1/6);  
    set(hAx(2),'Xlim',[1 size(dwi_vol,4)],'box','off','xtick',[],'ytick',get(hAx(2),'YLim'),'FontUnits','normalized','fontsize',1/6);              
    title(axes_bval_intensity_plot_top,'Bval (L) vs Signal Intensity (R)');
            
    % Plot translation and rotation - optional
    if ~isempty(movement_params_file) && movement_params_file.exist()
        % Load - only get for diffusion weighted images
        movement_params = movement_params_file.dlmread();
        movement_params = movement_params(bval_file.dlmread() ~= 0,:);
        
        % Do translation first
        trans_params = movement_params(:,1:3);
        plot(axes_plot_trans,trans_params(:,1),'LineWidth',1,'color','blue');
        hold(axes_plot_trans,'on');    
        plot(axes_plot_trans,trans_params(:,2),'LineWidth',1,'color','cyan');
        plot(axes_plot_trans,trans_params(:,3),'LineWidth',1,'color','magenta');
        set(axes_plot_trans,'Xlim',[1 size(trans_params,1)],'Ylim',[min(trans_params(:)) max(trans_params(:))],'box','off','YAxisLocation','right','xtick',[],'ytick',[min(trans_params(:)) max(trans_params(:))],'XColor',[0 0 0],'YColor',[0 0 0],'FontUnits','normalized','fontsize',1/6);              
        title(axes_plot_trans,'Translation');
        
        % Rotation
        rot_params = movement_params(:,4:6);
        plot(axes_plot_rot,rot_params(:,1),'LineWidth',1,'color','blue');
        hold(axes_plot_rot,'on');    
        plot(axes_plot_rot,rot_params(:,2),'LineWidth',1,'color','cyan');
        plot(axes_plot_rot,rot_params(:,3),'LineWidth',1,'color','magenta');
        set(axes_plot_rot,'Xlim',[1 size(rot_params,1)],'Ylim',[min(rot_params(:)) max(rot_params(:))],'box','off','YAxisLocation','right','xtick',[],'ytick',[min(rot_params(:)) max(rot_params(:))],'XColor',[0 0 0],'YColor',[0 0 0],'FontUnits','normalized','fontsize',1/6);       
        title(axes_plot_rot,'Rotation');
                
        % Add legend to rotation graph - note axes are not converted to RAS, 
        % so just leave the legend as "e1", "e2", and "e3"
        l = legend(axes_plot_rot,'e1','e2','e3','Orientation','Horizontal');
        set(l,'Position',[pos_plot_rot(1) pos_plot_rot(2)+pos_plot_rot(4)+3*padding pos_plot_rot(3) 2*padding]);
        drawnow
    end
       
    % Create stats text box   
    stats_string = {'','    Shells: '};
    bvals_unique = unique(bvals);
    for i = 1:length(bvals_unique)
        stats_string = horzcat(stats_string, ['        ' num2str(bvals_unique(i)) ': ' num2str(length(find(bvals == bvals_unique(i))))]); %#ok<AGROW>
    end    
    stats_string = horzcat(stats_string,{'','    Estimated sigma: ', ['        ' num2str(sigma_est)]});
    stats_height = pos_info(2)-(pos_footer(2)+pos_footer(4)+4*padding+chi_squared_area_height);
    stats_width = 0.65*(1 - (pos_chi_squared_plot(1)+pos_chi_squared_plot(3)) - padding);
    pos_stats_text = [1-stats_width-2*padding pos_info(2)-stats_height-2*padding stats_width stats_height];
    uicontrol('style','text','units','normalized','String',stats_string,...
        'FontUnits','Normalized','FontSize',1/35,...
        'Position',pos_stats_text,'HorizontalAlignment','left',...
        'BackgroundColor',[0.95 0.95 0.95],'Parent',f);
    
    % Create plots of masks
    pos_mask_text = [pos_mask_axial(1) pos_mask_axial(2)+pos_mask_axial(4)+padding mask_width mask_header];
    uicontrol('style','text','units','normalized','String',{'Masks'},...
        'FontUnits','Normalized','FontSize',1,'FontWeight','bold',...
        'Position',pos_mask_text,...
        'BackgroundColor',[0.85 0.85 0.85],'Parent',f);
    
    % get centroid of mask volume in RAS configuration    
    mask_vol_RAS = nifti_utils.load_untouch_nii_vol_scaled_RAS(mask_file.get_path(),'logical');  
    rp_mask_vol = regionprops(double(mask_vol_RAS),'Centroid','Area');
    centroid = round(rp_mask_vol.Centroid);
    if length(centroid) == 2 % vol might be 1D or 2D by matlabs standards if trailing dimensions are 1
        centroid(3) = 1;
    end
    centroid = centroid([2 1 3]); % (x,y,z) => (j,i,k)
    
    % Get visualizer    
    b0_vol_min = prctile(b0_vol(mask_vol),5);
    b0_vol_max = prctile(b0_vol(mask_vol),99);
    b0_vol = (b0_vol-b0_vol_min)./(b0_vol_max-b0_vol_min);
    dv = dwmri_visualizer({mask_vol,stats_mask_vol}, ...
                          b0_vol, ...
                          mask_vol, ...
                          xform_RAS, ...
                          'outlines', ...
                          {'r','g'; ...
                           0.25,0.25});
    % Plot axial
    dv.plot_slice(centroid(3),'axial','slice',[],axes_mask_axial);

    % Plot coronal 
    dv.plot_slice(centroid(2),'coronal','slice',[],axes_mask_coronal);

    % Plot sagittal
    dv.plot_slice(centroid(1),'sagittal','slice',[],axes_mask_sagittal);

    % set all to axis image
    axis(axes_mask_axial,'image');
    axis(axes_mask_coronal,'image');
    axis(axes_mask_sagittal,'image');
    drawnow
end
