function dti_stats_pdf_file = stats_plot(job_dir,dwmri_file,bval_file,mask_file,stats_mask_file,sigma_est,dt_file,v1_file,fa_file,outliers_file,chi_squared_plot,dwi_model_file,movement_params_file)

    disp('---');
    disp('Creating PDF...');
    
    % Create "PDF" directory
    pdf_dir = system_utils.directory(job_dir,'PDF');
    pdf_dir.mkdir();
       
    % Page 1: Plot RESTORE DTI statistics
    f_page1 = lib.page1_plot(dwmri_file, ...
                             bval_file, ...
                             mask_file, ...
                             stats_mask_file, ...
                             sigma_est, ...
                             outliers_file, ...
                             chi_squared_plot, ...
                             dwi_model_file, ...
                             movement_params_file);
        
    % Page 2: Plot DT
    f_page2 = lib.page2_plot(dwmri_file, ...
                             bval_file, ...
                             mask_file, ...
                             dt_file, ...
                             v1_file, ...
                             fa_file);
    
    % Save and merge PDFs
    dti_stats_pdf_file = system_utils.file(pdf_dir,'dti_stats.pdf');
    page1_pdf_file = system_utils.file(pdf_dir,'page1.pdf');
    page2_pdf_file = system_utils.file(pdf_dir,'page2.pdf');
    
    % Sometimes figures are saved before everything is "set". This might
    % fix it?
    drawnow
    
    print(f_page1,'-painters','-dpdf','-r600',page1_pdf_file.get_path());
    print(f_page2,'-opengl','-dpdf','-r600',page2_pdf_file.get_path());
    
    % Merge PDFs
    system_utils.system_with_errorcheck(['gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=' dti_stats_pdf_file.get_path() ' ' page1_pdf_file.get_path() ' ' page2_pdf_file.get_path() ],'Failed to merge output PDFs');

    % Remove single-paged pdfs
    page1_pdf_file.rm();
    page2_pdf_file.rm();
end
