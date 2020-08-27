function S = genorder(n,ordering,positive)
%GENORDER Generate colexicographical ordering
%   C = GENCOLEX(N) produces a sequence of (multi-)indices, starting from
%   (0,...,0) and up to (N1,...,Nd), sorted in the colexicographic order.
%
%   C = GENCOLEX(N,'sym') produces indices ranging from (-N1,...,-Nd) to
%   (N1,...,Nd).


d = numel(n);

switch d
    case 1
        S = (-n:n)';
    case 2
        [k1,k2] = ndgrid(-n(1):n(1),-n(2):n(2));
        switch ordering
            case 'colex'
                S = [k1(:),k2(:)];
            case 'lex'
                k1 = k1'; k2 = k2';
                S = [k1(:),k2(:)];
            case 'gcolex'
                S = [k1(:),k2(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
            case 'glex'
                k1 = k1'; k2 = k2';
                S = [k1(:),k2(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
        end
    case 4
        [k1,k2,k3,k4] = ndgrid(-n(1):n(1),-n(2):n(2),-n(3):n(3),-n(4):n(4));
        switch ordering
            case 'colex'
                S = [k1(:),k2(:),k3(:),k4(:)];
            case 'lex'
                k1 = permute(k1,[2 3 4 1]);
                k2 = permute(k2,[4 3 2 1]);
                k3 = permute(k3,[4 3 2 1]); % ou [1 2]?
                k4 = permute(k4,[4 3 2 1]);
                S = [k1(:),k2(:),k3(:),k4(:)];
            case 'gcolex'
                S = [k1(:),k2(:),k3(:),k4(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
            case 'glex'
                k1 = permute(k1,[2 3 4 1]);
                k2 = permute(k2,[4 3 2 1]);
                k3 = permute(k3,[4 3 2 1]); % ou [1 2]?
                k4 = permute(k4,[4 3 2 1]);
                S = [k1(:),k2(:),k3(:),k4(:)];
                [~,i] = sort(sum(abs(S),2));
                S = S(i,:);
        end
end

if positive
    S(~all(S>=0,2),:) = [];
end
            


end

