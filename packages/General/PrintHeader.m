function PrintHeader(proctime,opts)
    %Author: Feng-Yi Liao & Yang Zheng
    %        SOC Lab @UC San Diego    
    % Time setup and display
    proctime = toc(proctime);
    if opts.verbose
        [header,myline1,myline2] = Header();
        method = opts.method; % Set method to display
        fprintf(myline1)
        fprintf('Spectral Bundle Method for Primal SDPs - v1.0\n')
        fprintf('Authors: F. Liao, Y. Zheng @SOC Lab \n')
        fprintf(myline1)
        fprintf('Initilization is done in %.4f seconds.      \n',proctime);
        fprintf('Algorithm              : %s\n',method);

        fprintf('Free variables         : %i                \n',opts.K.f);
        fprintf('Non-negative variables : %i                \n',opts.K.l);
        fprintf('Second-order cones     : %i (max. size: %i)\n',length(find(opts.K.q ~=0)),max(opts.K.q));
        fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(find(opts.K.s ~=0)),max(opts.K.s));
        fprintf('Affine constraints     : %i                \n',opts.m);
        
        fprintf('Current eigenvectors   : %i                  \n',opts.current);
        fprintf('Penalty parameter      : %.2f                \n',opts.rho);
        fprintf('Past eigenvectors      : %i                \n',opts.past);
        
        fprintf(myline2);
        fprintf(header);
        fprintf(myline2);
    end
end