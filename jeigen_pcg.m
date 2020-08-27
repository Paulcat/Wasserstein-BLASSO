function [ B, trajectories, B_list ] = jeigen_pcg(As, varargin)
%# [ B trajectories ] = jeigen_pcg(As, varargin)
%# Given an m * m * Q set of matrices A_1, ..., A_Q  of size m * m,
%# find a m * m matrix B such that the matrices  B * A_q * inv(B)
%# are as diagonal as possible (sum of squares of off-diag terms )
%# - This works in the real or in the complex case
%# - It uses a precondition non linear conjugate gradient method
%#   (Polak-Ribiere style)

%# to do: find proper scalings for tol and regul

tstart = cputime;

[ m jk Q ] = size(As);
assert(m==jk, 'square matrices, please');

if Q==1 %# for Q=1, jeigen = eig
    [ B D ] = eig(As);
    B = inv(B);
    trajectories = [];
    return
end

%# parameters
defaults = {
    'maxiter', 500, ...
    'tol', 1e-6, ...
    'regul', 1e-2, ...
    'init', 'eye' , ...
    'verbose', true
    };
[ maxiter, tol, regul, init, verbose ] = ...
    processoptions (varargin, defaults{:} );

%# extra param
regul_on = 1e-6;
inner_balancing = false;

if verbose
    fprintf(1,'joint diagonlization of %i matrices of size %i\n', Q, m);
end

%# init: if a square matrix is passed as argument, use that.
if all(size(init)==[ m m ])
    B = init;
else
    switch init
        case 'eye'
            B = eye(m);
        case 'rand'
            B = randn(m) + 1i * randn(m);
        case 'algebra'
            [ B D ] = eig( As(:,:,1) - As(:,:,2) ); %# why not? because...
            B = inv(B);
        otherwise
            error('Unknown init option');
    end
end
As = updateA(As, B, inv(B));
o = off(As);

B_list = {};
B_list{end+1} = B;

if verbose
    fprintf(1,'JD crit. after init : %13.10f\n', o);
end

%# balance
[ L As iterbal ] = jd_balance_QN(As, false);
B = diag(L)*B;
o = off(As);

if verbose
    fprintf(1,'JD crit. after balancing powers : %13.10f\n', o);
end


%# main loop ################################################################
%# keep trajectories
jdcrit = zeros(maxiter,1);
moveoff  = zeros(maxiter,1);
moveon  = zeros(maxiter,1);
numbt  = zeros(maxiter,1);
gradoff  = zeros(maxiter,1);
gradon   = zeros(maxiter,1);
%#

G = zeros(m);
H = zeros(m);
D = zeros(m);
Hc = zeros(m);
T = zeros(m);
S = zeros(m);
R = zeros(m);

%# for CG reset
cg_reset = m*m; %
k = 0;

%# Init
[ G H Hc ] = der12(As); %# derivatives
R = -G;
S = applyprecon(R, H, Hc, regul, regul_on);
D = S;
delta_new = real(R(:)'*D(:));
g0 = 2 * real(G(:)'*D(:));  %# gradient in search direction

for iter=1:maxiter
    
    %# Line search ##
    [ As onew T nbt report ] = linesearch(As, o, D, g0) ;
    if (nbt==-1)
        if verbose, fprintf(1,'Max backtrack reached.\n'); end
        break
    end
    
    %# update
    B = T * B;
    gain = o - onew;
    o = onew;
    
    B_list{end+1} = B;
    
    [ G H Hc ] = der12(As); %# new derivatives
    gradoff(iter) = norm(G-diag(diag(G)), 'fro');
    gradon(iter)  = norm(diag(diag(G)  ), 'fro');
    
    %# stopping based on max relative move
    maxmove = max(abs(T-diag(diag(T))));
    if  ( maxmove < tol )
        if verbose, fprintf(1,'Tolerance reached: no relative move larger than %g\n', tol); end
        break
    end
    
    %# preconditionned conjugation with Polak-Ribiere
    R = -G;
    delta_old = delta_new;
    delta_mid = real(R(:)'*S(:));
    S         = applyprecon(R, H, Hc, regul, regul_on);
    delta_new = real(R(:)'*S(:));
    beta      = (delta_new - delta_mid) / delta_old;
    D         = S + beta * D;
    
    g0 = 2 * real(D(:)'* G(:));
    
    CG_RESET = (g0>0) | (beta<0) | (k == cg_reset); %#
    if CG_RESET
        D = S;
        k = 0;
        g0 = 2 * real(D(:)'* G(:));
        fprintf(1,'CG reset\n');
    else
        k = k + 1;
    end
    
    if inner_balancing  %# This step will interfer with the CG rationale...   Also changes g0 ?  How badly ?  Should we forget it?
        [ L As iterbal ] = jd_balance_QN(As, false);
        B = diag(L)*B;
    end
    
    %# if 0 # resetlarge
    %#   D /= max(abs(D(:)));
    %#   if verbose
    %# 	printf('Trim Large\n');
    %#   endif
    %# endif
    
    %# keeping track
    jdcrit(iter)  = onew;
    numbt(iter)   = nbt;
    moveon(iter)  = sqrt(sum( (abs(diag(T))-1).^2)); %# does not include the (optional) effect of balancing
    o2m_tmp_1 = (T-diag(diag(T)));
    moveoff(iter) = sqrt(sum( abs( o2m_tmp_1(:) ).^2));
    
    if verbose
        fprintf(1,'%4i|jd %18.13f (- %.1e)|on/off= %5.1e/%5.1e|grad on/off = %6.2e / %6.2e|nbt=%2i %s |CG: G0 = % 4.1e beta=% 4.1e |%g\n', ...
        iter, o, gain, moveon(iter), moveoff(iter), gradon(iter), gradoff(iter), nbt, report, g0, beta, CG_RESET  );
    end
    
end

if (iter==maxiter)
    if verbose
        fprintf(1,'Max number of iterations reached\n');
    end
end

if verbose
    dur =   cputime - tstart;
    fprintf(1,'JD criterion = % 25.22f.  Relative criterion = %.1e [in %i secs]\n', o, off_relative(As), ceil(dur));
end

%# final normalisation for one global scale and phase
%#B /= trace(abs(B));  ## Does not normalize the phase, but who cares?
B = B / (trace(B)); %# not normalize the phase, but who cares?
%# any better idea ?


%# Trim trajectories
keep = 1:(iter-1);
trajectories.moveon    = moveon  (keep);
trajectories.moveoff   = moveoff (keep);
trajectories.jdcrit    = jdcrit  (keep);
trajectories.numbt     = numbt   (keep);
trajectories.gradoff   = gradoff (keep);
trajectories.gradon    = gradon  (keep);



%###############################################################
function o = off(As)
[m jk Q] = size(As);
o = 0;
for q=1:Q
    M = As(:,:,q);
    M = M - (diag(diag(M)) );
    o = o + (M(:)'*M(:) );
end

%###############################################################
function or = off_relative(As)
%# return a normalized version of the JD criterion
[m jk Q] = size(As);
D = zeros(m,1);
M = zeros(m);
o = 0;
d = 0;
for q=1:Q
    M = As(:,:,q);
    D = diag(M);
    M = M - (diag(D) );
    o = o + ((M(:)'*M(:)) );
    d = d + (D'*D );
end
or = sqrt( o / d / m / Q );

%###############################################################
function A = updateA(A, T, iT)
Q = size(A,3);
for q=1:Q
    A(:,:,q) = T*A(:,:,q)*iT;
end

%###############################################################
function [ G H Hc ] = der12(AA)
%# 1st and 2nd-order quantities

[ m jk Q ] = size(AA);
G = zeros(m);
H = zeros(m);
Hc = zeros(m);
for q=1:Q
    A  = AA(:,:,q);
    dA = diag(A);
    
    Ao = A-diag(diag(A));
    G = G + (Ao*A' - A'*Ao );
    
    if nargout>1
        D  = (A .* conj(A))';
        r = sum(D,1);
        c = sum(D,2);
        H = H + (r + c - 2 * (D-diag(diag(D)) + real(dA*dA') ) ); %# broadcasting here
        Hc = Hc + (diag(r) + diag(c) - D - D' );
    end
    
end


%###############################################################
function S = applyprecon(R, H, Hc, regul, regul_on)

if 0 %# check rudimentary precon
    S = R / max(abs(H(:))) ;
    return
end

m = size(R,1);
S = R ./ (H + regul );
if 1
    e = (Hc + ones(m)/m + regul_on * eye(m) ) \ diag(R);
else
    e = diag(R) ./ diag(Hc + regul_on);
end
S = S + (diag( e / 2 - diag(S) ) );

%###############################################################
function [ Anew onew T ibt report ] =linesearch(As, o, D, g0)

%# line search seems robust wrt those parameters.
backtrackfactor  = 2;
maxbacktrack = 20;
c_wolfe  = 0.1; %

D2 = D*D;
m = size(D,1);
Anew = zeros(size(As));

%# 1: backtracking
alpha = 1;
btsuccess = false;
for ibt=1:maxbacktrack
    T = eye(m)  + alpha * D + 0.5 * (alpha*alpha)*D2; %#
    Anew = updateA(As, T, inv(T));
    onew = off(Anew);
    goodstep =  onew < (o + alpha *  c_wolfe * g0);
    if goodstep
        btsuccess = true;
        ibt = ibt - 1;
        break;
    end
    alpha = alpha / backtrackfactor;
end

if ~btsuccess
    ibt = -1;  %# used as a flag
    report = 'Backtrack failed';
    Anew = As;
    onew = o;
    T = eye(m);
    return;
else
    gain_bt = o-onew;
end

%# 2: then try one secant step, starting at backtrack value
Gnew = der12(Anew);
g1 = 2 * real(Gnew(:)'*D(:));
wrtbk = g0 / (g0-g1); %# step wrt backtrack step by the secant rule
alpha_sec = alpha * wrtbk;
Tsec = eye(m)  + alpha_sec * D + 0.5 * alpha_sec^2 *D2; %#
Asec = updateA(As, Tsec, inv(Tsec));
osec = off(Asec);
gainsec = o-osec;
%# check that secant step actually improves ovr backtrack (much needed during first iterations)
if gainsec > gain_bt  %# use secant result
    gain = gainsec;
    alpha = alpha_sec;
    Anew = Asec;
    onew = osec;
    T = Tsec;
    Gnew = der12(Anew);
    secimprove = gainsec / gain_bt;
else %# fall back to backtrack values
    gain = gain_bt;
    wrtbk = 1; %#
    secimprove = 0;
end
%# 3: should we iterate over secant steps ? In our experiments, one single step works really well.

%# relative gradient decrease during line search
%# this is computed `for information', so could be skipped (for a small relative speed).
%# It is not used as a stopping criterion for the line search (which is not iterative here)
g2 = 2 * real(   Gnew(:)'*D(:));
wc2 = g2 / g0;

report = sprintf('|LS: nbt=%2i g=%6.3f rstep=%5.3f wc2=% .1e', ...
    ibt, secimprove, wrtbk, wc2  );


%###############################################################
function [ L AA iter ] = jd_balance_QN(AA, verbose)
%# minimize the JD criterion wrt **diagonal** transforms
%# In other words, it fixes scales, which amounts to balancing
%# upper and lower energy.

if (nargin==1)
    verbose = false;
end

[ m jk Q ] = size(AA);

%# to deal with the global scale degeneracy
N = ones(m)/m;

tol = 1e-8;
itermax = 10;

%# this could be optimized loopless.  Worth it ?
B = zeros(m);
for q=1:Q
    Aq = AA(:,:,q);
    B = B + (Aq .* conj(Aq) );
end

L = ones(m, 1);
for iter=1:itermax
    
    o1 = sum(sum(B-diag(diag(B))));
    
    r = sum(B,1);
    c = sum(B,2);
    g = c(:)-r(:);
    H = diag(r) + diag(c) - B - B';
    d = - (H + N ) \ g;
    % *** bad HACK (by paul)!!!! *** %
    d(isnan(d)) = 0;
    %norm(-(H+N)*d - g,'fro')
    % **** %
    d = exp(d/2);
    L = L .* d;
    d = d .* d; %# square to act on B
    B = d .* B .*  (1.0 ./ d');  %# no need to update AA, only B is needed.
    
    ng = norm(g);
    if verbose
        fprintf(1,'Iter %i crit = %.6f ng = %.2e\n', iter, o1, ng);
    end
    
    if (ng < tol * o1)
        break
    end
    
end
%# now, we can update AA
for q=1:Q
    AA(:,:,q) = L .* AA(:,:,q) ./ L';
end


%###############################################################
function [varargout] =  processoptions  (params, varargin)
%# A simplified version of the parseparams code in octave by JFC
%# Original authors:
%# Author: Alexander Barth <abarth93@users.sourceforge.net>
%# Author: Aida Alvera Azcarate <aida@netecho.info>

%# process defaults
names     = varargin(1:2:end);
defvalues = varargin(2:2:end);
if (length(names) ~= length(defvalues))
    error ('needs defaults in the form of pairs keyword-value');
end

%# process specified values
pnames = params(1:2:end);
values = params(2:2:end);
if (length(pnames)~= length(values)) | ~iscellstr(pnames)
    error_as_caller ('options must be given as name-value pairs');
end

%# do the matching
varargout = defvalues;
isdefault = true(length(varargout), 1);
for ik=1:length(pnames)
    pname = pnames{ik};
    idx = find(strcmp(names, pname));
    if length(idx)==0
        error('unknown keyword %s', pname);
    end
    varargout{idx} = values{ik};
    isdefault(idx) = false;
end

v = find(strcmp(names, 'verbose'));
if (length(v)==1) && varargout{v}
    for ik=1:length(names)
        val =  varargout{ik};
        if isdefault(ik), defstr = 'default'; else defstr = ''; end
        if ischar(val)
            fprintf(1,'%14s = %14s [%s]\n', names{ik}, val, defstr);
        else
            if isscalar(val)
                fprintf(1,'%14s = %14g [%s]\n', names{ik}, val, defstr);
            else
                fprintf(1,'%14s = %14s [%s]\n', names{ik}, '<not a scalar>', defstr);
            end
        end
    end
end



