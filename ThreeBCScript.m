%% Considering Symmetry at r = 0 and Convection at r = R
% for i = 1:N
%     if i ==1
%         Main(i) = 1;
%         Upper(i) = -1;
%     elseif i == N
%         Lower(i) = -1;
%         Main(i) = 1 + htc*dr/K(i);
%     else
%         Lower(i) = -ALPHA(i)*dt/(dr^2) + ALPHA(i)*dt/(r(i)*dr);
%         Main(i) = 1 + 2*ALPHA(i)*dt/(dr^2) + rho_b*cp_b*W(i)*dt/(RHO(i)*CP(i));
%         Upper(i) = -ALPHA(i)*dt/(dr^2) - ALPHA(i)*dt/(r(i)*dr);
%     end
% end
% for n = 2:TS
%     Force(1) = 0;
%     Force(2:N-1)  = T(2:N-1,n-1) + rho_b*cp_b.*W(2:N-1)*dt*T_b./(RHO(2:N-1).*CP(2:N-1)) + QM(2:N-1)*dt./(RHO(2:N-1).*CP(2:N-1)) + q_mnp(2:N-1)*dt./(RHO(2:N-1).*CP(2:N-1));
%     Force(N) = htc*dr*T_inf/K(N);
%     T(:,n) = tridiag(Lower, Main, Upper, Force);
% end
% %%% Considering Symmetry at r = 0 and fixed Temperature at r = R
for i = 1:N
    if i ==1
        Main(i) = 1 + 2*ALPHA(i)*dt/(dr^2) + rho_b*cp_b*W(i)*dt/(RHO(i)*CP(i));
        Upper(i) = -2*ALPHA(i)*dt/(dr^2);
    elseif i == N
        Lower(i) = -0;
        Main(i) = 1;
    else
        Lower(i) = -ALPHA(i)*dt/(dr^2) + ALPHA(i)*dt/(r(i)*dr);
        Main(i) = 1 + 2*ALPHA(i)*dt/(dr^2) + rho_b*cp_b*W(i)*dt/(RHO(i)*CP(i));
        Upper(i) = -ALPHA(i)*dt/(dr^2) - ALPHA(i)*dt/(r(i)*dr);
    end
end

for n = 2:TS
    Force = zeros(N,1);
    Force(1:N-1)  = T(1:N-1,n-1) + rho_b*cp_b.*W(1:N-1)*dt*T_b./(RHO(1:N-1).*CP(1:N-1)) + QM(1:N-1)*dt./(RHO(1:N-1).*CP(1:N-1)) + q_mnp(1:N-1)*dt./(RHO(1:N-1).*CP(1:N-1));
    Force(N) = T_b;
    T(:,n) = tridiag(Lower, Main, Upper, Force);
end
% %%% Considering Symmetry at r = 0 and Adiabatic Boundary Condition at r = R
% for i = 1:N
%     if i ==1
%         Main(i) = 1;
%         Upper(i) = -1;
%     elseif i == N
%         Lower(i) = -1;
%         Main(i) = 1;
%     else
%         Lower(i) = -ALPHA(i)*dt/(dr^2) + ALPHA(i)*dt/(r(i)*dr);
%         Main(i) = 1 + 2*ALPHA(i)*dt/(dr^2) + rho_b*cp_b*W(i)*dt/(RHO(i)*CP(i));
%         Upper(i) = -ALPHA(i)*dt/(dr^2) - ALPHA(i)*dt/(r(i)*dr);
%     end
% end
%
% for n = 2:TS
%     Force = zeros(N,1);
%     Force(1) = 0;
%     Force(2:N-1)  = T(2:N-1,n-1) + rho_b*cp_b.*W(2:N-1)*dt*T_b./(RHO(2:N-1).*CP(2:N-1)) + QM(2:N-1)*dt./(RHO(2:N-1).*CP(2:N-1)) + q_mnp(2:N-1)*dt./(RHO(2:N-1).*CP(2:N-1));
%     Force(N) = 0;
%     T(:,n) = tridiag(Lower, Main, Upper, Force);
% end