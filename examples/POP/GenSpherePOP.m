%Generate POP with a unit sphere constraint
clc; clear;

for num = 1
    for n = [20]
        function_number = num;
        for z = 1:length(n)
            N = n(z);
            fprintf('num = %d,N = %d',num,N);
            %mpvar('x',[N,1]);%Remeber to use "mpvar". It is faster than "syms"
            x = sdpvar(N,1);
            function_name = "";
            f = 0;
            if function_number == 1 % "Broyden"
                f = ((3-2*x(1))*x(1)-2*x(2)+1)^2;
                for i =2:N-1
                    f = f + ((3-2*x(i))*x(i)-x(i-1)-2*x(i+1)+1)^2;
                end
                f = f+ ((3-2*x(N))*x(N)-x(N-1)+1)^2;
                f = f + sum(x).^2;
                function_name = "Broyden";
            elseif function_number ==2 % "Rosenbrock"
                f = 1;
                for i = 2:N
                    f = f + 100*(x(i)+x(i-1)^2)^2+(1-x(i))^2;
                end
                f = f + sum(x).^2;
                function_name = "Rosenbrock";
            elseif function_number ==3 % "Quartic"
                for i = 1:N
                    for j = i:N
                        f = f + x(i)*x(j) + x(i)^2*x(j) - x(j)^3 - x(i)^2*x(j)^2;
                    end
                end
                %f = f + sum(x).^2;
                function_name = "Quartic";
            elseif function_number ==4 % "random" 
                rng('default');
%                 c = rand(nchoosek(n+4,n),1);
%                 xx = monomials(x,[0,1,2,3,4]);
%                 f = c.'*xx;
                c = randn(nchoosek(n+4,n),1);
                v = monolist(x,4);
                f = c.'*v;
                function_name = "Random";
            else
                warning('The input function is not supported!');
            end
            
            R = 1; %radius of the shpere
            h = 0;
            for i = 1:N
                h = h + x(i)^2;
            end
            h = h - R;
            degree = 5;
            option.name = function_name + '_sphere' + num2str(N) + '_R' + num2str(R);
            gamma = sdpvar(1,1);
            %[s1,c1] = polynomial(x,2);
            phi_coe = sdpvar(nchoosek(n+2,n),1);
            s       = monolist(x,2);
            phi     = phi_coe.'*s;
            F = sos(f-gamma-phi*h);
            %solvesos(F,-gamma)
            [model,recoverymodel] = export(F,-gamma,sdpsettings('solver','sedumi'),[gamma;phi_coe]); %export the problem data in sedumi form
            At = model.A.';
            b  = model.b;
            c  = model.C;
            K  = model.K;
            save(option.name,'At','b','c','K');
            %findbound_modify(f,[],h,degree,option);
            %findbound(f,[],h,degree);
        end
    end
end