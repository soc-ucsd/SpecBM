function PrintHeader(proctime,opts)
% Time setup and display
    proctime = toc(proctime);
    if opts.verbose
        % Set method to display
        %method = 'SBMP';
        [header,myline1,myline2] = Header();
        method = opts.method;
%         if strcmpi(opts.solver,'hsde')
%             method = 'homogeneous self-dual embedding';
%         else
%             method = opts.solver;
%         end
        fprintf(myline1)
        fprintf('Spectral Bundle Method for Primal SDPs - v1.0\n')
        fprintf('Authors: F. Liao, Y. Zheng @SOC Lab \n')
        fprintf(myline1)
        fprintf('Initilization is done in %.4f seconds.      \n',proctime);
        fprintf('Algorithm              : %s\n',method);

%         if any(strcmpi(opts.solver,{'primal','dual'}))
%         fprintf('Adaptive penalty       : %i\n',opts.adaptive);
%         end
%         fprintf('Scale data             : %i\n',opts.rescale);
         fprintf('Free variables         : %i                \n',opts.K.f);
         fprintf('Non-negative variables : %i                \n',opts.K.l);
         fprintf('Second-order cones     : %i (max. size: %i)\n',length(find(opts.K.q ~=0)),max(opts.K.q));
         fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(find(opts.K.s ~=0)),max(opts.K.s));
         fprintf('Affine constraints     : %i                \n',opts.m);

         fprintf('Penalty parameter      : %.2f                \n',opts.rho);
         fprintf('Current eigenvectors   : %i                  \n',opts.current);
         fprintf('Past eigenvectors      : %i                \n',opts.past);

         fprintf(myline2);
         fprintf(header);
         fprintf(myline2);
    end
end