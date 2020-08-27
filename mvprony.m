function [supp,amp,info] = mvprony(mm,n,d,options)
%MVPRONY Multivariate Prony extraction
%   [S,A] = MVPRONY(MM,N,D,options) extracts the support S and amplitudes A
%   of a discrete measure explaining the moment matrix MM. N specifies the
%   maximal order contained in MM and D the dimension.
%
%   Options:
%       'factorized' - whether MM is factorized [ (0) | 1 ]
%       'shift_mode' - expression for multiplication matrices [ harmouch |
%           klep | (kunis) ]
%       'jdiag_step' - procedure for joint diagonalization [ (cardoso) |
%           random ]
%       'ordering' - monomials ordering [ lex | (colex) | glex | gcolex ]
%       'tol' - tolerance for svd (1e-3)

M = size(mm,1);

% set options
factorized = getoptions(options,'factorized',0);
shift_mode = getoptions(options,'shift_mode','kunis');
jdiag_step = getoptions(options,'jdiag','cardoso');
ordering   = getoptions(options,'ordering','colex');
tol        = getoptions(options,'tol',1e-3);

% only 'positive' moments?
positive = (M==prod(n+1));

% generate ordering
o  = genorder(n,ordering,positive);
o1 = genorder(n-1,ordering,positive);

% extract submatrix
[~,i1] = ismember(o1,o,'rows');

% compute shifting operators
Shift = cell(1,d);
for i=1:d
    oi = o1;
    oi(:,i) = oi(:,i)+1; % shift index along i-th dimension
    
    [~,ishift] = ismember(oi,o,'rows'); % position of shifted index in matrix (indexed by o)
    Shift{i} = sparse(i1,ishift,ones(length(i1),1),M,M);
end

% *** compute multiplication matrices ***
% ***************************************
N = cell(1,d);
switch shift_mode
    case 'kunis'
        if factorized
            [U,S,V] = mysvdf(mm(i1,:),tol);
            
            for i=1:d
                mi   = Shift{i}*mm; mi = mi(i1,:);
                N{i} = U' * mi * mm(i1,:)' * V * diag(1./S);
            end
        else
            [U,S,V] = mysvd(mm(i1,i1),tol);
            
            for i=1:d
                Mi   = Shift{i}*mm; Mi = Mi(i1,i1);
                N{i} = U' * Mi * V * diag(1./S);
            end
        end
        
    case 'klep'
        if factorized
            [U,S,V] = mysvdf(mm(i1,:),tol);
            
            for i=1:d
                mi   = Shift{i}*mm; mi = mi(i1,:);
                N{i} = diag(sqrt(1./S)) * U' * mi * mm(i1,:)' * V * diag(sqrt(1./S)); 
            end
        else
            [U,S,V] = mysvd(mm(i1,:),tol);
            
            for i=1:d
                Mi   = Shift{i}*mm; Mi = Mi(i1,i1);
                N{i} = diag(sqrt(1./S)) * U' * Mi * V * diag(sqrt(1./S));
            end
        end
        
    case 'harmouch'
        if factorized
            [U,S,V] = mysvdf(mm(i1,:),tol);
            
            for i=1:d
                mi   = Shift{i}*mm; mi = mi(i1,:);
                N{i} = diag(1./S) * U' * mi * mm(i1,:)' * V;
            end
        else
            [U,S,V] = mysvd(mm(i1,i1),tol);
            
            for i=1:d
                Mi   = Shift{i}*mm; Mi = Mi(i1,i1);
                N{i} = diag(1./S) * U' * Mi * V;
            end
        end    
end
s = length(S); % number of reconstructed points

% in theory, multiplication matrices should commute if hierarchy collapses
if d==2
    e = norm(N{1}*N{2}-N{2}*N{1}, 'fro') / norm(N{1}*N{2}, 'fro');
    fprintf('Relative commutation error: %.3f\n', e);
end

% *** joint diagonalization ***
% *****************************
if strcmp(jdiag_step,'random')
    lambda = getoptions(options,'lambda',rand(1,d));
    lambda = reshape(lambda,[1 1 d]);
    
    Nco = sum(lambda .* cat(3, N{:}), 3);
    [H,h] = eig(Nco); h = diag(h);
else
    As = cell2mat(N);
    As = reshape(As,[size(N{1}),d]);
    [H,~] = jeigen_pcg(As,'init','eye');
    H = inv(H);
end

% *** compute support (from eigenvalues of multiplication matrices) ***
% *********************************************************************
supp    = zeros(s,d);
modulus = zeros(s,d);
for i=1:d
    ei = diag(inv(H)*N{i}*H);
    supp(:,i) = mod(-angle(ei)/(2*pi),1);
    
    % mean modulus may isolate outliers
    modulus(:,i) = abs(ei);
end
modulus = mean(modulus,2);

% *** compute amplitudes (Vandermonde system) ***
% ***********************************************
% matching index range to moment range: [0 n]->[-n n], [-n n]->[-2n 2n]
if positive
    ofull = genorder(n,ordering,0);
else
    ofull = genorder(2*n,ordering,1);
end
% estimated matrix
Gr    = reshape(ofull,[size(ofull,1),1,d]); % frequencies sorted following ordering
Sr    = reshape(supp,[1,s,d]);
Fsupp = exp(-2i*pi* (sum(Gr.*Sr,3)) ); % broadcasting

[~,ids] = marginals(n,ordering,positive); % all moments (no repetition)
if factorized
    [i,j] = ind2sub(M,ids);
    c     = mm(i,:) .* conj(mm(j,:));
    c     = sum(c,2);
else
    c = mm(ids);
end

% amplitudes: solve least-square system
amp = real(Fsupp \ c(:));

% *** additional infos ***
% ************************
info.modulus = modulus;
end

function [U,S,V] = mysvd(M,tol)
disp('computing svd...')

[U,S,V] = svd(M,'econ'); S = diag(S);
I = find(abs(S)>tol*max(abs(S)));
U = U(:,I); V = V(:,I); S = S(I);

disp('done')
end

function [U,S,V] = mysvdf(X,tol)
disp('computing svd...')

[U,S] = svd(X,'econ'); S = diag(S).^2;
I = find(abs(S)>tol*max(abs(S)));
U = U(:,I); V = U; S = S(I);

disp('done')
end
