function [P, Pc, mask, n_s] = set_decomp (n, m, n_o, as_type)
% Gratton Zhang 2017
% n number of variables
% m number of subspaces
% n_o size of the overlap
% as_type 'as' or 'ras' or ?? for weighting

P = [];
Pc = [];
mask = [];
n_s = zeros(m,1) + NaN;

if (n < m)
    error('The number of subspaces is too large.');
end
%if (n < (m+1)*n_o )
if (n < m*n_o )
    error('The overlap is too large.');
end

% Define the decomposition. Each subspace has n_s variables,
% beginning with x(n_beg) and endding with x(n_end). 
r = (n-n_o) - m*(floor((n-n_o)/m));
n_s(1:m) = floor((n-n_o)/m) + n_o;
n_s(floor((m-r)/2)+1 : floor((m-r)/2)+r) = n_s(floor((m-r)/2)+1 : floor((m-r)/2)+r) + 1; 

n_s_cum = [0; cumsum(n_s(1:m-1))];
n_beg = n_s_cum(1:m) - n_o.*(0:m-1)' + 1;
n_end = n_beg + n_s - 1; 

II = speye(n,n);
shared_2 = [];

for i = 1:m 
    mask{i} = ones(n_s(i), 1);
    if (strcmp('as', as_type))
    elseif (strcmp('ras', as_type))
        if (i >= 2)
            mask{i}(1 : floor(n_o/2)) = 0;
        end
        if (i <= m-1)
            mask{i}(n_s(i)-(n_o-floor(n_o/2))+1 : n_s(i)) = 0;
        end
    else
        if (i >= 2)
            mask{i}(1 : n_o) = 0.5;
        end
        if (i <= m-1)
            mask{i}(n_s(i)-n_o+1 : n_s(i)) = 0.5;
        end
    end
    P{i} = II(:,n_beg(i) : n_end(i));
    if (i <= m-1)
        shared_2 = [shared_2; (n_beg(i+1):n_end(i))'];
    end
end

if (n_o >= 2) % Serge: please check whether the definition of Pc is correct or not.
    l_shared_2 = length(shared_2);
    %Pc = interp1(shared_2,speye(l_shared_2,l_shared_2),1:n,'linear',0);
    %Pc = interp1([0;shared_2;n+1],[zeros(1,l_shared_2);speye(l_shared_2,l_shared_2);zeros(1,l_shared_2)],1:n,'linear',0);
    %Pc = interp1(shared_2,eye(l_shared_2,l_shared_2),1:n,'linear',0);
    %keyboard
    %Pc = interp1([0;shared_2;n+1],[zeros(1,l_shared_2);eye(l_shared_2,l_shared_2);zeros(1,l_shared_2)],1:n,'linear',0);
end




