function y = pathhomology(D,pmax)

% NB. For a more computationally efficient and capable implementation
% building on this one, see another function by the same name (and its
% helpers) that takes vargin, was written by M. Yutin, and is (or was
% originally) packaged nearby. The relative benefits of the present
% implementation are simplicity and a permissive license.
%
% Constructs the path complex of a digraph D up to dimension pmax, and its
% concomitant (reduced) path homology up to dimension pmax-1.
% For background, details, etc. see
%     [PPHDN] Chowdhury, S. and Memoli, F. "Persistent path homology of
%     directed networks." SODA (2018).
%     SIAM version:	https://doi.org/10.1137/1.9781611975031.75
%     ACM mirror:	https://dl.acm.org/citation.cfm?id=3175345
%     arXiv:        https://arxiv.org/abs/1701.00565
% or https://arxiv.org/abs/1207.2834 and other papers on the subject by
% Grigor'yan, Yong, Muranov, Yau, and collaborators. 
%
% The reduced and non-reduced Betti numbers differ only in dimension 0. The
% non-reduced zeroth Betti number is just the number of weakly connected
% components of D, which is easily obtained via the command
% conncomp(D,'Type','weak'). We thank Samir Chowdhury for this observation
% (see also 3.6 of https://arxiv.org/abs/1207.2834).
%
% The use of non-regular versus regular path homology is motivated by
% example 3.14 of https://arxiv.org/abs/1207.2834 in light of potential
% applications where 2-cycles seem worthy of "hole" status, but note that
% our preference is precisely opposite to the paper's.
%
% Since ranks are computed using svds, this code can give wrong answers at
% times. If in doubt, trust, but verify using comments in the code that can
% be uncommented to use symbolic matrices and rank (if the symbolic toolbox
% is available). Note that this is much more computationally expensive.
%
% The output is a struct with fields ordered by dimension(+1). Fields
% include
%     allPaths:     all(owed) paths of a given length (measured by edges)
%     omega:        matrix representation of invariant paths
%     bdry:         boundary maps of the underlying chain complex
%     betti:        Betti numbers of the chain complex
%
% Copyright (c) 2020, BAE Systems. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
% contributors may be used to endorse or promote products derived from this
% software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Preliminaries
n = size(D.Nodes,1);
% If D.Nodes.Name does not exist, make a simple lex ordered one
if ~any(cellfun(@numel,strfind(fieldnames(D.Nodes),'Name')))
    D.Nodes.Name = cellstr(dec2base((1:n)',10,ceil(log10(n))));
end
% Enforce lex order for convenience
D = reordernodes(D,sort(D.Nodes.Name));
% Warn if there are any loops and then remove them
loops = findedge(D,1:n,1:n);
if any(loops)
    warning('removing loops!');
    D = rmedge(D,1:n,1:n);
end
% Set any weights to unity
if any(cellfun(@numel,strfind(fieldnames(D.Edges),'Weight')))
    D.Edges.Weight = ones(size(D.Edges,1),1);
end
% Adjacency and identity matrices
A = adjacency(D);

%% Compute number of paths of length <= pmax between each pair of vertices
numPaths = speye(n);    % paths of length 0
numLen_p = speye(n);
for p = 1:pmax
    numLen_p = numLen_p*A;
    numPaths = numPaths+numLen_p;
end

%% Construct all(owed) paths of length <= pmax between pairs of vertices
% This could be done merely by taking powers of a symbolic version of the
% adjacency matrix, but that would require a lot of heavy lifting under the
% hood and not be amenable to low-level implementations. 
allPaths = cell(1,pmax+1);    % double entendre
[S,T] = find(numPaths);
for i1 = 1:numel(S)
    s = S(i1);
    t = T(i1);
    % paths_st is a cell containing paths of length <= pmax from s to t
    paths_st = reashortestpaths(D,s,t,numPaths(s,t));
    pp1 = cellfun(@numel,paths_st); % pp1 = p+1;
    % Organize paths by length
    for i2 = 1:numel(pp1)
        allPaths{pp1(i2)} = [allPaths{pp1(i2)};paths_st{i2}];
    end
end

%% Lex order paths of a given length
for p = 0:pmax
    allPaths{p+1} = sortrows(allPaths{p+1});
end

%% Construct the chain complex
indA = cell(1,pmax+1);  % indices of allowed p-paths in elementaries
bdryA = cell(1,pmax+1); % modified boundary map
omega = cell(1,pmax+1); % space of invariant p-paths
bdry = cell(1,pmax+1);  % invariant boundary map
for p = 0:pmax
    if numel(allPaths{p+1})
        %% Indices of allowed p-paths in set of elementary p-paths
        % Everything is lex ordered, use fancy radix manipulation below:
        indA{p+1} = 1+(allPaths{p+1}-1)*n.^(p:-1:0)';
        % NB. Analytically determining indices of regular p-paths seems
        % hard. While we can apparently do this in the ambient set of lex
        % ordered elementary p-paths, it is intricate (see commented code
        % at the end of this file), and currently based on occult
        % numerology (i.e., we have no proof that the technique works,
        % though presumably such a proof must be inductive). Yet unanswered
        % is the question of how to index allowed p-paths in the ambient
        % set of regular p-paths.
        %% Initialize modified boundary map
        % bdryA sends allowed p-paths to elementary (in particular,
        % regular) p-1 paths: we do more lex ordering cuteness here 
        bdryA{p+1} = sparse(n^p,size(allPaths{p+1},1));
        %% Populate columns of the modified boundary map
        for i1 = 1:size(allPaths{p+1},1)
            curPath = allPaths{p+1}(i1,:);
            for i2 = 0:p
                % Possibly faster alternative to kill i2th vertex than
                % temp = curPath(setdiff(1:(p+1),i2+1));
                temp = curPath([1:i2,(i2+2):(p+1)]);
                % More fancy radix manipulation here makes this look
                % scarier than it really is
                bdryA{p+1}(1+(temp-1)*n.^((p-1):-1:0)',i1) = (-1)^i2;
            end
        end
        %% Remove rows corresponding to allowed (p-1)-paths and get kernel
        % This yields the the space of invariant p-paths
        if p > 0
            bdryTemp = bdryA{p+1}(setdiff(1:n^p,indA{p}),:);
        else
            bdryTemp = sparse(1,n);
        end
        % Kill zero columns before computing kernel to save LOTS of time
        omega{p+1} = null(full(bdryTemp(any(bdryTemp,2),:)));
        % % Uncomment to use symbolic/exact calculation
        % omega{p+1} = null(sym(bdryTemp(any(bdryTemp,2),:)));
        %% Construct invariant boundary maps
        % % Uncomment to use symbolic/exact calculation
        % if ~numel(omega{p+1})
        %     omega{p+1} = sym(zeros(size(bdryA{p+1},2),0)); 
        % end
        temp = bdryA{p+1}*omega{p+1};
        if p > 0
            bdry{p+1} = temp(indA{p},:);
        else
            bdry{p+1} = temp;
        end
    else
        bdryA{p+1} = sparse(n^p,0);
        omega{p+1} = 0;
        bdry{p+1} = sparse(n^p,0);
    end
end

%% Compute dimensions of boundary map images and kernels
dim_im = zeros(1,pmax+1);
dim_ker = zeros(1,pmax+1);
for p = 0:pmax
    % Kill zero columns before computing rank
    dim_im(p+1) = rank(full(bdry{p+1}(any(bdry{p+1},2),:)));
    % dim_im(p+1) = rank(sym(bdry{p+1}(any(bdry{p+1},2),:)));
    dim_ker(p+1) = size(bdry{p+1},2)-dim_im(p+1);   % rank-nullity theorem
end

%% Compute Betti numbers
betti = dim_ker(1:pmax)-dim_im(2:end);

%% Output
y.allPaths = allPaths;
y.bdryA = bdryA;
y.omega = omega;
y.bdry = bdry;
y.betti = betti;

%% Indices of regular p-paths in set of elementary p-paths (lex ordered)
% % n = 4;
% % pmax = 6;
% %% Indices of regular p-paths in set of elementary p-paths (lex ordered)
% % A touch of fancy radix stuff lets one get indices of allowed p-paths in
% % the set of elementary p-paths. Thus the missing piece is how to get
% % indices of allowed p-paths in the set of regular p-paths...
% indD = repmat([ones(1,n-1),2],[1,n-1])  % index differences
% indD = indD(1:end-1)
% indDs = cell(1,n-1) % will hold suffixes of indD
% Xi = repmat([1,2],[1,n-1])	% entries are intercalated, hence Greek Xi
% Xi = Xi(1:end-1)
% % Indices of regular p-paths in elementary p-paths, starting with p = 2.
% % For p = 0 there is a single index (1) and for p = 1 the indices are 1:n
% indR = [];
% indR{1} = 1:n
% indR{2} = setdiff((1:(n^2)),1:(n+1):(n^2))
% %% Loop to get
% for p = 3:pmax
%     %% Get suffixes of indD
%     len = ((n-1):-1:1)*(n-1)^(p-2)-1;
%     for j = 1:(n-1)
%         indDs{j} = indD((end-len(j)+1):end);
%     end
%     %% Get intercalating differences
%     Xi = n*Xi+(-1).^((1:(2*n-3))+p)-(n-2);
%     %% Assemble next indD
%     indD = [];
%     for j = 1:(n-1)
%         if mod(j,2)
%             indD = [indD,indDs{1},Xi(j)];
%         else
%             indD = [indD,fliplr(indDs{1}),Xi(j)];
%         end
%         indDs = fliplr(indDs(2:end));
%     end
%     indD = [indD,fliplr(indD(1:end-1))];
%     %% Get initial offset and recover index from differences
%     indR{p} = cumsum([indR{p-1}((n-1)^(p-2)+1),indD]);
%     %% Brute force for comparison and verifying assertion
%     elem = dec2base(0:(n^p-1),n,p);
%     ind = find(all(diff(elem,1,2),2))';
%     reg = elem(ind,:);
%     assert(all(indR{p}==ind),'foobar');
% end