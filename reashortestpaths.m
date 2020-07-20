function p = reashortestpaths(D,s,t,K)

% Given a weighted digraph D, source vertex s, target vertex t, and
% positive integer K, this returns the K shortest paths in D from s to t.
%
% The recursive enumeration algorithm (REA) of Jimenez and Marzal
% (https://doi.org/10.1007/3-540-48318-7_4) is used instead of the
% algorithms of Eppstein or Aljazzar and Leue since it performs well in
% practice, "it relies on quite simple data structures, and it can be
% easily implemented". We borrow much but not all of their description and
% notation in lieu of documenting this code more extensively.
%
% NB. There is a minor error in the REA paper: in Step B.1, it is necessary
% to require that $\pi^1(u)$ actually exists, or in other words that that
% $u$ is actually reachable from $s$.
%
% The implementation deliberately uses cells for simplicity vs elegance.
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

%% Notation
n = size(D.Nodes,1);
m = size(D.Edges.EndNodes,1);
p = cell(1,n);  % p{v}{k} is the kth shortest path from s to v
C = cell(1,n);  % C{v} is a SET of candidate next-shortest paths

%% If weights do not already exist, initialize them to unity
if ~any(strcmp(D.Edges.Properties.VariableNames,'Weight'))
    D.Edges.Weight = ones(m,1);
end

%% Step A.1
% Compute p{v}{1} for all v...and set k = 1
k = 1;
tree = shortestpathtree(D,s);
for v = 1:n
    p{v}{k} = shortestpath(tree,s,v);
    C{v} = {};
end

%% Step A.2
% Repeat until p{t}{k} does not exist or no more paths are needed (or in
% other words, while p{t}{k} exists and more paths are needed):
while ~isempty(p{t}{k}) && k < K
    %% Step A.2.1
    % Set k = k+1 and compute p{t}{k} (and update C) by calling
    % nextPath(D,s,t,p,C)
    k = k+1;
    [p,C] = nextPath(D,s,t,p,C);
end

%% Output
if isempty(p{t}{end})
    warning('there are not K paths in D from s to t');
end
p = p{t};

end
%% End main function

%% Local function
function [p,C] = nextPath(D,s,v,p,C)

k = numel(p{v})+1;

%% Step B.1
if k == 2
    % Initialize a SET of candidates to the next shortest path from s to v.
    % As a proxy for a SET, we use a cell whose entries are unique by
    % construction (cf. Step B.5). 
    pre = predecessors(D,v);
    for i = 1:numel(pre)
        u = pre(i);
        if ~isempty(p{u}{1})    % fixes minor error in the REA paper
            temp = [p{u}{1},v];
            if ~isequal(temp,p{v}{1})
                C{v}{end+1} = temp;
            end
        end
    end
end

%% Steps B.2-B.5
if v == s && k == 2
    %% Step B.2
    % Go to Step B.6
else
    %% Step B.3
    % Let u and kk be the node and index, respectively, such that p{v}{k-1}
    % = [p{u}{kk},v]
    u = p{v}{k-1}(end-1);
    cloneB3 = repmat({p{v}{k-1}(1:end-1)},size(p{u},1),size(p{u},2));
    kk = find(cellfun(@isequal,p{u},cloneB3));
    % Could also do the above with the following elementary approach:     
    %     for kk = 1:numel(p{u})
    %         if isequal(p{v}{k-1},[p{u}{kk},v])
    %             break;  
    %         end
    %     end

    %% Step B.4
    % If p{u}{kk+1} has not already been computed, then compute it by
    % calling nextPath(D,s,u,p,C)
    if numel(p{u}) < kk+1
        [p,C] = nextPath(D,s,u,p,C);
    end
    
    %% Step B.5
    % If p{u}{kk+1} exists, then insert [p{u}{kk+1},v] in the SET C{v}
    if ~isempty(p{u}{kk+1})
        cloneB5 = repmat({[p{u}{kk+1},v]},size(C{v},1),size(C{v},2));
        if ~any(cellfun(@isequal,C{v},cloneB5))
            % [p{u}{kk+1},v] is not already in C{v}, so append it
            C{v}{end+1} = [p{u}{kk+1},v];
        end
        % Could also do the above with the following elementary approach:    
        %     append = 1;
        %     for i = 1:numel(C{v})
        %         if isequal(C{v}{i},[p{u}{kk+1},v])
        %             append = 0;
        %             break;
        %         end
        %     end
        %     if append
        %         C{v}{end+1} = [p{u}{kk+1},v];
        %     end
    end
end

%% Step B.6
% If ~isempty(C{v}), then select and delete the path with minimum length
% from C{v} and assign it to p{v}{k}, else p{v}{k} does not exist
if ~isempty(C{v})
    %% Compute lengths of paths in C{v}
    len = zeros(1,numel(C{v}));
    for i = 1:numel(C{v})
        edges = findedge(D,C{v}{i}(1:(end-1)),C{v}{i}(2:end));
        len(i) = sum(D.Edges.Weight(edges));
    end
    %% Select the path with minimum length from C{v}
    % Note that we are free to tiebreak in this way per discussion in the
    % REA paper since we use a shortest path tree in Step A.1
    ind = find(len==min(len),1,'first');
    %% Delete this path from C{v} and assign it to p{v}{k}  
    p{v}{k} = C{v}{ind};
    C{v} = C{v}(setdiff(1:numel(C{v}),ind));
else
    p{v}{k} = [];	% does not exist
end

end