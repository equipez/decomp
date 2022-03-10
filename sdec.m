function [x, fx, exitflag, output] = sdec(fun, x0, options)

verbose=0;

n = length(x0);

% Default options:
dim = 1; % Dimension of decomposition
m = 2; % Number of subspaces
n_o = 2;
as = 'ras';
cs = 'rdfs';
% Coarse spaces: 
% '': no coarse space
% 'pss': previous subspace steps (not good)
% 'rpss': restricted previous subspace steps (OK)
% 'pfs': previous full space step (good)
% 'dfs': decomposition of previous full space step (not good)
% 'rdfs': restricted decomposition of previous full space step (good,
% similat to 'pfs')
% 'cg': conjugate gradient
% 'grad': current gradient (OK, not very good)
% 'shem': Spectral Harmonically Enriched Multiscale coarse space (only
% for minmal surface/Laplace and 2D decomposition)
ls = false;
% For 'pss', 'rpss', 'pfs', 'dfs', 'rdfs', 'cg', whether take Last Successful
% step (ls=true) or simply last step (ls=false). No significant
% difference observed. ls = fale even a bit better in a few tests. Strange!!! 
maxit = 1000;
ftarget = -Inf;
debug = true;
print = 0;

Delta0 = 1;

%tol_Delta = 1e-6; tol_Pred = 1e-12; tol_d = 1e-8; tol_cg = 1e-8;
tol_g = 1e-6; tol_Delta = 0; tol_Pred = 0; tol_d = 0; tol_cg = 1e-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 3)
    if (~isa(options, 'struct') && ~isempty(options))
        error('sdec:InvalidOptions', 'options should be a structure or empty.');
    end
    if (isfield(options,'dim')) 
        dim = max(options.dim, 1);
    end
    if (isfield(options,'m')) 
        m = max(options.m, 1);
    end
    if (isfield(options,'n_o')) 
        n_o = max(options.n_o, 0);
    end
    if (isfield(options,'as')) 
        as = options.as;
    end
    if (isfield(options,'ls')) 
        ls = options.ls;
    end
    if (isfield(options,'cs')) 
        cs = options.cs;
    end
    if (isfield(options,'Delta0'))
        Delta0 = options.Delta0;
    end
    if (isfield(options,'tol_g')) 
        tol_g = max(options.tol_g, eps);
        tol_cg = max(eps, min(1e-2*tol_g, 1e-10));
    end
    if (isfield(options,'tol_Delta')) 
        tol_Delta = max(options.tol_Delta, eps);
    end
    if (isfield(options,'maxit')) 
        maxit = int32(options.maxit);
    end
    if (isfield(options,'ftarget')) 
        ftarget = options.ftarget;
    end
    if (isfield(options,'print')) 
        print = options.print;
    end
    if (isfield(options,'debug')) 
        debug = options.debug;
    end
end


if (~isa(fun, 'function_handle') && ~isa(fun, 'char'))
    error('sdec:InvalidFun', 'fun should be a function handle or function name.');
end

if (~isnumeric(x0) || ~isvector(x0))
    error('sdec:InvalidX0', 'x0 should be a numerical vector or scalar.')
end
x0 = x0(:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work with function handels instread of function names to avoid using 'feval'.
if (isa(fun, 'char'))
    fun = str2func(fun);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim == 2
    m1 = m(1);
    m2 = m(2);
    m = m1*m2;
end

if (m <= 1)
    P{1} = speye(n,n);
    mask{1} = ones(n,1);
    Pc = [];
    cs = '';
end
              
if (dim == 1 && m > 1)
    [P, ~, mask] = set_decomp(n, m, n_o, as);
elseif (dim == 2 && m > 1)
    if strcmpi(cs, 'shem') 
        coarse = 8; % Strange code inherited from old version; affects only set_decomp2
%        [P, Pc, mask, ratio_ff] = set_decomp2(sqrt(n), sqrt(m), n_o, as, coarse, options.stride);
        [P, Pc, mask, ratio_ff] = set_decomp2(sqrt(n), [m1,m2], n_o, as, coarse, options.stride);
        %Qc = qr(Pc, 0);
        Qc = Pc; % No qr seems better 
    else
        coarse = 0;
%        [P, Pc, mask, ratio_ff] = set_decomp2(sqrt(n), sqrt(m)0, n_o, as, coarse, 1);
        [P, Pc, mask, ratio_ff] = set_decomp2(sqrt(n), [m1,m2], n_o, as, coarse, 1);
        Qc = [];
    end
end
Predc = 0;

eta0 = 0;
rholimit = 1;
if ((strcmpi(as, 'ras') || strcmpi(as, 'wras')) && n_o > 0) 
    if (dim == 1 && m > 1) || (dim == 2 && m1*m2 > 1 && min(m1, m2) == 1)
        rholimit = 1/2;
    elseif dim == 2 && m1 > 1 && m2 > 1
        rholimit = 1/4; 
    end
end

% Possibilities:
%eta1 = 0.1*rholimit; eta2 = 0.9*rholimit;
%eta1 = 0.05*rholimit; eta2 = rholimit;
%eta1 = 0.25*rholimit; eta2 = 0.75*rholimit;
%eta1 = 0.125*rholimit*0.75; eta2 = 0.375*rholimit*0.75;  
%eta1 = 0.125*rholimit; eta2 = 0.375*rholimit;
%eta1 = 0.1*rholimit; eta2 = 0.75*rholimit;
%eta1 = 0.1*rholimit; eta2 = 0.5*rholimit; %%%
%eta1 = 0.05*rholimit; eta2 = 0.5*rholimit;
%eta1 = 0.01*rholimit; eta2 = 0.5*rholimit;
%eta1 = 0.1*rholimit; eta2 = rholimit; %%%
%eta1 = 0.1*rholimit; eta2 = 0.9*rholimit; %%%
eta1 = 0.1*rholimit; eta2 = 0.75*rholimit; %%%
%eta1 = 1.1*rholimit; eta2 = inf; %%% volating the theory

% Possibilities:
%gamma1 = 0.5; gamma2 = 4;
%gamma1 = 0.25; gamma2 = 2; %%%
gamma1 = 0.5; gamma2 = 2;

% Possibilities:
% Note that the w_k in the paper is 1-wcs
wcs = 0.25; % Generally good 
%wcs = 0.5; % For 'pfs', big wcs seems favorable
%wcs = 0.75;

ghist = NaN(maxit+1, 1);
fhist = NaN(maxit+1, 1);
Dhist = NaN(maxit+1, 1);
rhist = NaN(maxit+1, 1);

x = x0;
Delta = Delta0;
Dhist(1) = Delta0;

[fx, g, H] = fun(x);
f0 = fx;
ng0 = norm(g);
fhist(1) = fx;
ghist(1) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exitflag = -1;

if verbose
    disp('    rho         nd            delta          mdec         fdec          ||g||          f')
end

succ = 0;


for iter = 1:maxit
    d_sync = zeros(n,1);
    Pred =0; 

    if exist('Pc_old', 'var')
        switch cs
        case {'rpss', 'pss', 'rdfs', 'dfs', 'pfs'}
            Qc = Pc_old; % Note that Pc_old may be a matrices with orthogonal collumns, and hence orthogonalization will produce the same result as normalization 
        case 'cg'
            Qc = [Pc_old, -g];
        case 'grad'
            Qc = [-g];
        otherwise
            Qc = [];
        end
        normQ = sqrt(sum(Qc.^2, 1));
        indexQ = (normQ >= 0.01*max([normQ,tol_g]));
%        Qc = Qc(:,indexQ);
        Qc = Qc(:,indexQ)*diag(normQ(indexQ).^(-1));
    elseif ismember(cs, {'grad', 'cg'})
        Qc = [-g/norm(g)];
    elseif ~strcmpi(cs, 'shem')
        Qc = [];
    end
     %Qc = qr(Qc, 0); % orthogonalization  
     % Other Possibilities:
     %Qc = [-g];

    ng = 0;
    for i = 1:m
       ng = ng + (norm(P{i}'*g))^2; 
    end
    ng = sqrt(ng);

    nd = 0;
    for i = 1:m
        dl = zeros(length(P{i}(1, :)), 1);
        gl = P{i}'*g; 
        Hl = P{i}'*H*P{i};
        Deltal = (norm(gl)/ng)*Delta; 
        [dl, Predl, mess] = tcg (gl, Hl, Deltal, length(gl), tol_cg);

        if Predl ~= Predl % tcg may fail
            dl = zeros(size(gl));
            Predl = 0;
        end
        assert (Predl >= 0);

        nd = nd + dl'*dl;

        d_sync = d_sync + P{i}* (dl .* mask{i});

        Pred = Pred + Predl;

        if (strcmpi(cs, 'rpss'))
            Pc(:,i) = P{i} * (dl.*mask{i} ); 
        end
        if (strcmpi(cs, 'pss'))
            Pc(:,i) = P{i} * dl;
        end
    end
    if exist('Qc', 'var') && ~isempty(Qc)
        gc = Qc'*g;
        Hc = Qc'*H*Qc;
        %Deltac = 0.5*Delta;
        %Deltac = Delta/norm(full(Qc)); % norm(M) is not available for sparse M
        Deltac = Delta; 
        % Scaling Delta with norm(Qc) does not make much difference for
        % all the coarse spaces except SHEM, for which not scaling is better. 
        [dc, Predc, msg, k] = tcg (gc, Hc, Deltac,  length(gc), tol_cg);
        dc = Qc*dc;
        %wcs = max(0, min(0.99, d_sync'*(d_sync - dc)/((d_sync-dc)'*(d_sync - dc))))
        d_trail = (1-wcs)*d_sync + wcs*dc;
        
        Pred = (1-wcs)*Pred + wcs*Predc;
        %nd = nd + dc'*dc;  %??? % Theoretically not reasonable
    else
        d_trail = d_sync;
    end
    % Possibilities: 
    %nd = norm(d_trail); % Theoretically not reasonable
    nd = sqrt(nd);

    if (norm(d_trail) <= max(eps, tol_d*Delta0))
        exitflag = 3;
        break;
    end

    x_trail = x + d_trail;
    [f_trail,g_trail,H_trail]=fun(x_trail);

%    if (Pred <= max(eps, tol_Pred*(1+abs(f0))))
%        Pred
%        exitflag = 4;
%        break;
%    end
    Ared = fx-f_trail;
%    rho = (Ared+eps)/(Pred+eps);
    rho = Ared/Pred;
    if rho ~= rho
        exitflag = 4;
        break;
    end
    rhist(iter) = rho;
    
    if (~ls || rho > eta0 || succ == 0)
        switch cs
        case {'rpss', 'pss'}
            Pc_old = Pc;
        case 'rdfs'
            for i = 1:m
                Pc_old(:,i) = P{i}*(P{i}'*d_trail.*mask{i});
            end 
        case 'dfs'
            for i = 1:m
                Pc_old(:,i) = P{i}*(P{i}'*d_trail);
            end 
        case {'cg', 'pfs'}
            Pc_old = d_trail;
        end
    end
    if (rho > eta0)
        succ = 1;
        x = x_trail;
        fx = f_trail;
        g = g_trail;
        H = H_trail;
    end
    fhist(iter + 1) = fx;
    ghist(iter + 1) = norm(g)/ng0;

    if (fx <= ftarget)
        exitflag = 0;
        break;
    end
    if (norm(g) <= max(eps, tol_g*ng0))
        exitflag = 1;
        break;
    end

    if (rho < eta1)
        % Possibilities:
        %Delta = gamma1*Delta; 
        Delta = gamma1*nd;
    %elseif (rho > eta2 && nd >= Delta/gamma2)
    elseif (rho > eta2)% && Delta <= 1000*norm(g))
        % Possibilities:
        %Delta = gamma2*Delta;
        %Delta = gamma2*nd; % This seems to be very important. 
        Delta = max(Delta, gamma2*nd);
    else
        %Delta = nd;
        Delta = Delta;
    end
    Dhist(iter + 1) = Delta;

    if (Delta <= max(eps, tol_Delta*Delta0))
        exitflag = 2;
        break;
    end

if verbose
%  disp(sprintf('%d %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e', iter, rho, delta, value, f_old-f_trial, norm(g),f_val) );
%   disp(sprintf('%d %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e', iter, rho, nd, Delta, Pred, fx, norm(g),fx, Pred, Predc) );
   disp(sprintf('%d %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e %+13.11e ', iter, rho, nd, Delta, Pred, Ared, norm(g),fx) );
end
end

output.fhist = fhist(1 : iter);
output.Dhist = Dhist(1 : iter);
output.rhist = rhist(1 : iter);
output.ghist = ghist(1 : iter);
output.iter = iter;
output.relg = norm(g)/ng0;

if (print == 1)
    disp(sprintf('   flag = %d - rel gradient = %+13.11e - iteration = %d', exitflag, output.relg, output.iter))
end

function [ s, value, msg, k, gplus, interior ] = tcg (g, H, delta,  kmax, epsi);
% Example:
% Gratton Zhang 2017
% H=diag([1 2 3]);
% g = [1;1;1]; delta=20; kmax=10; epsi=1e-6;
% [ s, msg, k, value, gplus, interior ] = tcg (g, H, delta,  kmax, epsi);
% delta = 1;
% [ s, msg, k, value, gplus, interior ] = tcg (g, H, delta,  kmax, epsi);
%
% Initializations
%
value = 0;

ng0 = norm(g);

%  Define TCG adapt

% example of termination rule  
 eps_term = epsi * ng0  ;
 eps_term = min( epsi, sqrt( ng0 ) ) * ng0  ;

s  = 0*g;
v  =  g ;
p  =  -v;
%
for k = 1:kmax
   Hp    = H * p;
   kappa = p' * Hp;

   % Check if negative curvature encountered

   if ( kappa <= 0 )
      sigma    = boundary_pt( s, p, delta);
      s        = s + sigma * p;
      msg      = [ 'st_tcg - negative curvature - ', int2str(k), ' iterations' ];
      value    = value - sigma * g' * p - 0.5 * sigma^2 * kappa;
      gplus    = g + sigma * Hp;
      interior = 0;
      return
   end

   gv    = g' * v;
   alpha = gv / kappa;

   % Check if outside the trust region
   if ( ( norm( s+alpha*p) ) >= delta )
      sigma    = boundary_pt( s, p, delta );
      s        = s + sigma * p;
      msg      = [ 'st_tcg - solution on the boundary - ', int2str(k), ' iterations '];
      value    = value - sigma * g' * p - 0.5 * sigma^2 * kappa;
      gplus    = g + sigma * Hp;
      interior = 0;
      return
   end
%
%  CG body
%

   s      = s + alpha * p;
   value  = value + 0.5 * alpha * gv; 
   gplus  = g + alpha * Hp;
   vplus  =  gplus ;
   gpvp   = gplus' * vplus;
   beta   = gpvp / gv;
   p      = -vplus + beta * p;
   g      = gplus;
   v      = vplus;

%  termination ?
            
   ngplus = sqrt( gpvp );
   if  ( ngplus < eps_term )
      msg=[ 'st_tcg - interior solution - ', int2str(k), ' iterations' ];
      interior = 1;
      return
   end
end
msg=['st_tcg - no solution in ',int2str(k),' iterations'];
interior = 1;

function sigma = boundary_pt( s, p, delta );
if delta ==0
  sigma = 0;
else
   np2   = p' * p;
   ns2  = s' * s;
   smp   = s' * p;
   sigma =(-smp+sqrt(smp^2+np2*(delta^2-ns2)))/(np2);
end
