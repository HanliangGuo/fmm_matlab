%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FMM structure for a toy problem (2D Laplace problem)
% algorithm follows Greengard & Rokhlin 1987
% regular (non-adaptive) quadtree
% particles are in a unit box centered at (0.5 ,0.5)
%
% Hanliang Guo
% Apr 28 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

% initialize
N = 1e4;    % number of source/target points
z = rand(N,1)+1i*rand(N,1);
q = ones(N,1);  % strength of charges;
p = 5;  % multipole terms
n = 3;  % level of refinement

VAL = 1;    % validation flag

tic;
base_ID_global = -1;
disp('Initialization');
for l = 1:n
    base_ID_global = base_ID_global+4^(l-1);
    for ibox = 1:4^l
        ibox_global = base_ID_global + ibox;
        Box(ibox_global).level = l;
        Box(ibox_global).ID = ibox;
        Box(ibox_global).center = getboxcenter(l,ibox);
        Box(ibox_global).pID = getparentID(l,ibox);
        Box(ibox_global).cID = getchildrenID(l,ibox);
        Box(ibox_global).nlist = getneighbors(l,ibox);
        Box(ibox_global).ilist = getinteractionlist(l,ibox);
        Box(ibox_global).a = zeros(1,p+1);
        Box(ibox_global).b = zeros(1,p+1);
        Box(ibox_global).c = zeros(1,p+1);
        Box(ibox_global).d = zeros(1,p+1);
        Box(ibox_global).particlelist = []; % this step may need to be performed at every time step
    end
end
% whos Box
%%
% Preparation
% Construct a matrix of useful binomial coefficients
nCk = zeros(2*p,p);
for i_n = 1:2*p+1
    for i_k = 1:min(i_n,p+1)
        nCk(i_n,i_k) = nchoosek(i_n-1,i_k-1);
    end
end

% Upward Pass
% Step 1
% Form multipole expansions of potential field due to particles in each box
% about the box center at the finest mesh level
disp('Step 1, S2M');
% tic;
idx_remaining = (1:N)';
for ibox_global = 1:4
    boxsize = 1/4;
    zC = Box(ibox_global).center;
    list1 = find(abs(real(z-zC))<boxsize);
    list2 = find(abs(imag(z-zC))<boxsize);
    idx = intersect(list1, list2);
    Box(ibox_global).particlelist = [Box(ibox_global).particlelist idx'];

end
for l = 2:n
    for ibox=1:4^l
        ibox_global = sum(4.^(1:l-1))+ibox;
        pID_global  = Box(ibox_global).pID+sum(4.^(1:l-2));
        particlecandidate = Box(pID_global).particlelist;
        boxsize = 1/2^(l+1);
        zC = Box(ibox_global).center;
        list1 = find(abs(real(z(particlecandidate)-zC))<boxsize);
        list2 = find(abs(imag(z(particlecandidate)-zC))<boxsize);
        idx = Box(pID_global).particlelist(intersect(list1,list2));
        Box(ibox_global).particlelist = idx;
        if l == n
%     % Form a p-term multipole expansion Phi(n,ibox), by using Theorem 2.1 (S2M)
%     % phi(z) = Q log(z)+\sum_{k=1}^p a_k/z^k
%     % Q = sum_{i=1}^m q_i; a_k = -sum_{i=1}^m q_i^z_i^k/k;
            Box(ibox_global).a(1) = Box(ibox_global).a(1) + sum(q(idx));
            for k = 1:p
                Box(ibox_global).a(k+1) = Box(ibox_global).a(k+1) - sum((q(idx).*(z(idx)-zC).^k))/k;
            end
            
        end
    end
end
% toc

% Step 2
% Form multipole expansions about the centers of all boxes at all coarser
% mesh levels, each expansion representing the potential field due to all
% particles contained in one box.
disp('Step 2, M2M');
% tic;
for l = n-1:-1:1
    for ibox=1:4^l
        % Form a p-term multipole expansion Phi(l,ibox) by using Lemma 2.3 to shift the
        % multipole expansion of each child's box center to the parent's
        % box center and adding them together
        ibox_global = sum(4.^(1:l-1))+ibox;
        zM = Box(ibox_global).center;
        
        cID = getchildrenID(l,ibox);
        ibox_global_c = sum(4.^(1:l))+cID;
        
        for ichild = 1:4
            Box(ibox_global).b(1) = Box(ibox_global).b(1) + Box(ibox_global_c(ichild)).a(1);
            zC = Box(ibox_global_c(ichild)).center;
            
            for ll = 1:p    
                Box(ibox_global).b(ll+1) = Box(ibox_global).b(ll+1) - Box(ibox_global_c(ichild)).a(1)*(zC-zM)^ll/ll;
                for k = 1:ll
                    Box(ibox_global).b(ll+1) = Box(ibox_global).b(ll+1) + Box(ibox_global_c(ichild)).a(k+1)*(zC-zM)^(ll-k)*nCk(ll,k);
                end
            end
        end
        Box(ibox_global).a = Box(ibox_global).b;
    end
end
% toc
%%
% Downward pass
% In the downward pass, interactions are consistently computed at the
% coarsest possible level. For a given box, this is accomplished by
% including interactions with those boxes which are well separated and
% whose interactions have not been accounted for at the parent's level.

% Step 3
% Form a local expansion about the center of each box at each mesh level
% l<=n-1. This local expansion describes the field due to all particles in
% the system that are not contained in the current box or its nearest
% neighbors. Once the local expansion is obtained for a given box, it is
% shifted, in the second inner loop to the centers of the box's children,
% forming the initial expansion for the boxes at the next level.

% Set Psihat(1,:) = 0;
disp('Step 3, Local expansion (M2L and L2L)');
% tic;
for l = 1:n-1
    for ibox = 1:4^l
        % Form Psi(l,ibox) by using lemma 2.4 to convert the multipole expansion Phi(l,j)
        % of each box j in interaction list of box ibox to a local
        % expansion about the center of box ibox, adding these local
        % expansions together, and adding the result to Psihat.
        % M2L
        ibox_global = sum(4.^(1:l-1))+ibox;
        Box(ibox_global).c = Box(ibox_global).d;
        iboxlist = Box(ibox_global).ilist;
        if isempty(iboxlist)==0     %interaction list is not empty
            zL = Box(ibox_global).center;
            for i = 1:length(iboxlist)
                ibox_global_m = sum(4.^(1:l-1))+iboxlist(i);
                zM = Box(ibox_global_m).center;
                Box(ibox_global).c(1) = Box(ibox_global).c(1) + Box(ibox_global_m).b(1)*log(zL-zM);
                for k = 1:p
                    Box(ibox_global).c(1) = Box(ibox_global).c(1) + Box(ibox_global_m).b(k+1)/(zM-zL)^k*(-1)^k;
                end
                
                kk = 1:p;
                pref = Box(ibox_global_m).b(kk+1)./(zM-zL).^kk.*(-1).^kk;
                for ll = 1:p
                    Box(ibox_global).c(ll+1) = Box(ibox_global).c(ll+1) + sum(pref.*nCk(ll+kk,ll+1)'./(zM-zL)^ll) - Box(ibox_global_m).b(1)/ll/(zM-zL)^ll;
                end
            end
        end
        
    end
    

    for ibox = 1:4^l
        % Form the expansion Psihat(l+1,j) for ibox's children by using lemma 2.5 to expand
        % Psi(l,ibox) about the children's box centers.
        % L2L
        ibox_global = sum(4.^(1:l-1))+ibox;
        zL = Box(ibox_global).center;
        cID = getchildrenID(l,ibox);
        ibox_global_c = sum(4.^(1:l))+cID;
        for i = 1:4
            ibox_global_ci = ibox_global_c(i);
            zT = Box(ibox_global_ci).center;
            for ll = 0:p
                for k = ll:p
                    Box(ibox_global_ci).d(ll+1) = Box(ibox_global_ci).d(ll+1) + Box(ibox_global).c(k+1)*nCk(k+1,ll+1)*(zT-zL)^(k-ll);
                end
            end
        end
    end
    
end
% toc;

% Step 4
% Compute interactions at finest mesh level
disp('Step 4, M2L at finest mesh level')
% tic;
for ibox = 1:4^n
    % Form Psi(l,ibox) by using lemma 2.4 converting the multipole expansion Psi(l,j) of
    % each box j in interaction list of box ibox to a local expansion about
    % the center of box ibox, adding these local expansions together, and
    % adding the result to Psihat(l,ibox).
    ibox_global = sum(4.^(1:n-1))+ibox;
    Box(ibox_global).c = Box(ibox_global).d;

    iboxlist = Box(ibox_global).ilist;
    if isempty(iboxlist)==0     %interaction list is not empty
        zL = Box(ibox_global).center;
        for i = 1:length(iboxlist)
            ibox_global_m = sum(4.^(1:n-1))+iboxlist(i);
            zM = Box(ibox_global_m).center;
            Box(ibox_global).c(1) = Box(ibox_global).c(1) + Box(ibox_global_m).a(1)*log(zL-zM);
            for k = 1:p
                Box(ibox_global).c(1) = Box(ibox_global).c(1) + Box(ibox_global_m).a(k+1)/(zM-zL)^k*(-1)^k;
            end
            
            kk = 1:p;
            pref = Box(ibox_global_m).a(kk+1)./(zM-zL).^(kk).*((-1).^kk);
            for ll = 1:p
%                 Box(ibox_global).c(ll+1) = Box(ibox_global).c(ll+1);
                
                temp = sum(pref./(zM-zL).^(ll).*nCk(ll+kk,ll+1)');
                Box(ibox_global).c(ll+1) = Box(ibox_global).c(ll+1) + (temp) - Box(ibox_global_m).a(1)/ll/(zM-zL)^ll;
            end
        end
    end
end
% toc;
% Local expansions at finest mesh level are now available. They can be used
% to generated the potential or force due to all particles outside the
% nearest neighbor boxes at finest mesh level.

% Step 5
% Evaluate local expansions at particles positions.
disp('Step 5, Evaluate local expansions at particles positions')
% tic;
Phi = zeros(N,1);
for ibox = 1:4^n
    % For every particle p_j located at the point z_j in box ibox, evaluate
    % Phi(n,ibox)(z_j)
    ibox_global = sum(4.^(1:n-1))+ibox;
    zC = Box(ibox_global).center;
%     For all particles in the box
    for i = Box(ibox_global).particlelist
%     local expansion about zC
        for ll = 0:p
            Phi(i) = Phi(i) + (Box(ibox_global).c(ll+1))*(z(i) - zC)^ll;
        end
    end
end

% Step 6
% Compute potential (or force) due to nearest neighbors directly
disp('Step 6, Evaluate potential due to nearest neighbors directly');
% tic;
Phi_direct0 = zeros(N,1);
Phi_direct1 = zeros(N,1);
for ibox = 1:4^n
    % For every particle p_j in box ibox, compute interactions with all
    % other particles within the box and its nearest neighbors.
    ibox_global = sum(4.^(1:n-1))+ibox;
    neighbors = Box(ibox_global).nlist;
    particles_inbox = Box(ibox_global).particlelist;
    particles_neighbors = [];
    for i = neighbors
        ibox_global_n = sum(4.^(1:n-1))+i;
        particles_neighbors = [particles_neighbors Box(ibox_global_n).particlelist];
    end
    
    temp = log(transpose(z(particles_inbox))-z(particles_inbox));
    diag_ind = 1+(length(particles_inbox)+1)*(0:(length(particles_inbox)-1));
    temp(diag_ind) = 0;
    Phi_direct0(particles_inbox) = sum(q(particles_inbox).*temp);
    Phi_direct1(particles_inbox) = sum(q(particles_neighbors).*log(transpose(z(particles_inbox))-z(particles_neighbors)));

end
Phi_direct = Phi_direct0 + Phi_direct1;
% toc;


% Step 7
disp('Step 7, Combine direct and far-field terms together');
% tic;
Potential = real(Phi_direct) + real(Phi);
elapsedTime = toc;
disp(['FMM time spent: ' num2str(elapsedTime) ' seconds']);

% Validation
if VAL == 1
    disp('Validation: Evaluate potential directly')
    tic;
%     Potential_test = zeros(N,1);
%     for i = 1:N
%         for j = 1:N
%             if i~=j
%                 Potential_test(i) = Potential_test(i) + q(j)*(log(z(i)-z(j)));
%             end
%         end
%     end

    Potential_test = q'.*log(transpose(z)-z);
    diag_ind0 = 1+(N+1)*(0:(N-1));
    Potential_test(diag_ind0) = 0;
    Potential_test = sum(Potential_test,2);

    Potential_test = real(Potential_test);
    % Potential_test = Potential_test - mean(Potential_test);
    elapsedTime = toc;
    disp(['Direct evaluation time spent: ' num2str(elapsedTime) ' seconds']);
    error1 = Potential_test-Potential;

    rms_error = norm(error1./Potential_test,2)/sqrt(N)
    rlt_error = norm(error1./Potential_test,2)/norm(Potential_test,2)
end

% Each local expansion is described by the coefficients of a p-term
% polynomial. Direct evaluation of this polynomial at a point yields the
% potential. But, the force is immediately obtained from the derivative
% which is available analytically. There is no need for numerical
% differentiation. Furthermore, due to the analyticity of Psi', there exist
% error bounds for the force of exactly the same form as (2.4), (2.11), and
% (2.15) in Greengard & Rokhlin 1987.

function pID = getparentID(l,ibox)
if ibox > 4^l
    disp('getparentID: boxID and level not consistent!');
    stop;
end
% pID = ceil(ibox/4^(l-1));
ibox_x = mod(ibox-1,2^l)+1;
ibox_y = ceil(ibox/2^l);
ibox_x_p = ceil(ibox_x/2);
ibox_y_p = ceil(ibox_y/2);
pID = (ibox_y_p-1)*2^(l-1) + ibox_x_p;
end

function cID = getchildrenID(l,ibox)
if ibox > 4^l
    disp('getchildrenID: boxID and level not consistent!');
    stop;
end
% cID = (ibox-1)*4 + (1:4)';
ibox_x = mod(ibox-1,2^l)+1;
ibox_y = ceil(ibox/2^l);
ibox_x1 = (ibox_x-1)*2+1;
ibox_x2 = (ibox_x-1)*2+2;
ibox_y1 = (ibox_y-1)*2+1;
ibox_y2 = (ibox_y-1)*2+2;
cID = [(ibox_y1-1)*2^(l+1) + ibox_x1;
    (ibox_y1-1)*2^(l+1) + ibox_x2;
    (ibox_y2-1)*2^(l+1) + ibox_x1;
    (ibox_y2-1)*2^(l+1) + ibox_x2;];
% showlist(l+1,[],cID);
end

function xC = getboxcenter(l,ibox)
if ibox > 4^l
    disp('getboxcenter: boxID and level not consistent!');
    stop;
end
ibox_x = mod(ibox-1,2^l)+1;
ibox_y = ceil(ibox/2^l);
xC = (ibox_x-.5)/2^l + 1i*(ibox_y-.5)/2^l;
end

function nlist = getneighbors(l,ibox)
nlist = [];
ibox_x = mod(ibox-1,2^l)+1;
ibox_y = ceil(ibox/2^l);

if ibox_y>1
    if ibox_x>1
        nID = (ibox_y-2)*2^l + ibox_x-1;    % bottom-left neighbor
        nlist = [nlist;nID];
    end

    nID = (ibox_y-2)*2^l + ibox_x;    % bottom neighbor
    nlist = [nlist;nID];
    if ibox_x<2^l
        nID = (ibox_y-2)*2^l + ibox_x+1;    % bottom-right neighbor
        nlist = [nlist;nID];
    end
end

if ibox_x>1
    nID = (ibox_y-1)*2^l + ibox_x-1;    % left neighbor
    nlist = [nlist;nID];
end

if ibox_x<2^l
    nID = (ibox_y-1)*2^l + ibox_x+1;    % right neighbor
    nlist = [nlist;nID];
end

if ibox_y<2^l
    if ibox_x>1
        nID = (ibox_y)*2^l + ibox_x-1;    % top-left neighbor
        nlist = [nlist;nID];
    end

    nID = (ibox_y)*2^l + ibox_x;    % top neighbor
    nlist = [nlist;nID];

    if ibox_x<2^l
        nID = (ibox_y)*2^l + ibox_x+1;    % top-right neighbor
        nlist = [nlist;nID];
    end
end
% showlist(l,ibox,nlist);
end

function ilist = getinteractionlist(l,ibox)
    ilist = [];
    pID = getparentID(l,ibox);
    p_nlist = getneighbors(l-1,pID);
    for i = 1:length(p_nlist)
        cID = getchildrenID(l-1,p_nlist(i));
        ilist = [ilist;cID];
    end
    nlist = getneighbors(l,ibox);
    ilist = setdiff(ilist,nlist);
%     showlist(l,ibox,ilist);
end

function showlist(l,ibox,list)
C = ones(2^l);
if ibox>0
ibox_x = mod(ibox-1,2^l)+1;
ibox_y = ceil(ibox/2^l);
C(ibox_y,ibox_x) = 30;
end
for i = 1:length(list)
    iibox = list(i);
    iibox_x = mod(iibox-1,2^l)+1;
    iibox_y = ceil(iibox/2^l);
    C(iibox_y,iibox_x) = 60;
end
image(C);
axis equal
end