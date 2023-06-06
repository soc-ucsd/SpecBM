function Out = SBM(At,b,c,K,userOpts)
%Main function
%Spectral bundle method Primal and Dual formulation 
%
%Syntax:
%
%
%Solve a standard semidefinite program in the following form
%
%        min   <C,X>               max  <b,y>
%  (P)   s.t.  A(X) = b       (D)  s.t. A^Ty+Z = C
%              X \in PSD                Z\in PSD
%                                      
%                                             

    if ~isfield(userOpts,'solver')
        userOpts.solver = "primal";% default
    end
  
    if userOpts.solver == "primal"
        Out = SBMP(At,b,c,K,userOpts);
    elseif userOpts.solver == "dual"
        Out = SBMD(At,b,c,K,userOpts);
    else
        fprintf('Incorrect solver input \n');
        return;
    end

end