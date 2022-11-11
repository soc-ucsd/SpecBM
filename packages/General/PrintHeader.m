function PrintHeader(proctime,opts)
% Time setup and display
    proctime = toc(proctime);
    if opts.verbose
        % Set method to display
        %method = 'SBMP';
        method = opts.method;
%         if strcmpi(opts.solver,'hsde')
%             method = 'homogeneous self-dual embedding';
%         else
%             method = opts.solver;
%         end
        fprintf('Initilization is done in %.4f seconds.      \n',proctime);
        fprintf('Algorithm              : %s\n',method);

%         if any(strcmpi(opts.solver,{'primal','dual'}))
%         fprintf('Adaptive penalty       : %i\n',opts.adaptive);
%         end
%         fprintf('Scale data             : %i\n',opts.rescale);
%         fprintf('Free variables         : %i                \n',K.f);
%         fprintf('Non-negative variables : %i                \n',K.l);
%         fprintf('Second-order cones     : %i (max. size: %i)\n',length(find(K.q ~=0)),max(K.q));
%         fprintf('Semidefinite cones     : %i (max. size: %i)\n',length(find(K.s ~=0)),max(K.s));
%         fprintf('Affine constraints     : %i                \n',opts.m);
%         if any(strcmpi(opts.solver,{'primal','dual','hsde'}))
%         fprintf('Consensus constraints  : %i                \n',sum(accumarray(Ech,1)));
%         else
%         fprintf('Nonorthogonal dimension: %i                \n',opts.sos.NonOrth);
%         end
        [header,myline1,myline2] = Header();
        fprintf(myline1);
        fprintf(header);
        fprintf(myline2);
    end
end