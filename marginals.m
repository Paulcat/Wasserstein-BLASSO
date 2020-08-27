function [ids,all_ids] = marginals(n,ordering,positive)
%MARGINALS Indices of marginal moments in moment matrix
%   IDS = MARGINALS(N,ORDER,POS) produces a cell containing the linear
%   indices in the moment matrix of moments for each marginal, without
%   repetition, sorted following ORDER
%
%   IDS,ALL_IDS = MARGINALS(N,ORDER,POS) also returns the indices of all
%   moments (ie 'diagonal' values), without repetition.
%
%   POS specifies whether moments range from -N to N (POS==0) or from 0 to
%   N (POS==1).

d = length(n);
o = genorder(n,ordering,positive);

% determine size of moment matrix
if positive
    M = prod(n+1);
else
    M = prod(2*n+1);
end

pows = reshape(o,[M,1,d]) - reshape(o,[1,M,d]);

% all moments indices in moment matrix (without repetition)
powlist = reshape(pows,M^2,d); % list of powers
disp('computing all indices');
[~,I] = unique(powlist,'rows','stable');
disp('done!')
all_moms = powlist(I,:); % all powers, without repetition
[Y,X] = meshgrid(1:M);
all_ids = [X(:),Y(:)]; % all indices, with repetition
all_ids = all_ids(I,:); % without repetition
all_ids = sub2ind([M M],all_ids(:,1),all_ids(:,2)); % linear indexing


% marginal moments indices in moment matrix
for i=1:d
    powsi = pows(:,:,i);
    [k,l] = find(powsi==0);
    
    I = k==1 | l==1; % avoid repetition
    i1 = k(I);
    i2 = l(I);
    
    ids{d-i+1} = sub2ind([M M],i1,i2); % linear indexing
end

% sort according to ordering
ofull = genorder(2*n,ordering,0);

% for all moments
[~,jall] = ismember(all_moms,ofull,'rows'); % position wrt ordering
[~,sort_all] = sort(jall);
all_ids = all_ids(sort_all,:);

% for moments of marginals
for i=1:d
        powsi = pows(:,:,i);
        moms  = zeros(length(ids{i}),d);
        moms(:,i) = powsi(ids{i}); % powers of i-th marginal moments
        
        [~,jmarg] = ismember(moms,ofull,'rows'); % position wrt ordering
        [~,sort_marg]  = sort(jmarg);
        ids{i} = ids{i}(sort_marg);
end





end

