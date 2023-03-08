function [ Aout ] = estimate_airlight( img, Amin, Amax, N, spacing, K, thres )

if ~exist('thres','var') || isempty(thres), thres = 0.01 ; end;
if ~exist('spacing','var') || isempty(spacing), spacing = 0.02 ; end; 
if ~exist('n_colors','var') || isempty(N), N = 1000 ; end; 
if ~exist('K','var') || isempty(K), K = 40 ; end; 

if ~exist('Amin','var') || isempty(Amin), Amin = [0,0.05,0.1]; end;
if ~exist('Amax','var') || isempty(Amax), Amax = 1; end;

if isscalar(Amin), Amin = repmat(Amin,1,3); end 
if isscalar(Amax), Amax = repmat(Amax,1,3); end

[img_ind, points] = rgb2ind(img, N);
[h,w,~] = size(img);

disp(img_ind);

idx_in_use = unique(img_ind(:));
disp(idx_in_use);
idx_to_remove = setdiff(0:(size(points,1)-1),idx_in_use);
points(idx_to_remove+1,:) = [];
img_ind_sequential = zeros(h,w);
for kk = 1:length(idx_in_use)
    img_ind_sequential(img_ind==idx_in_use(kk)) = kk;
end

[points_weight,~] = histcounts(img_ind_sequential(:),size(points,1));
points_weight = points_weight./(h*w);
if ~ismatrix(points), points = reshape(points,[],3); end 

angle_list = reshape(linspace(0, pi, K),[],1);

directions_all = [sin(angle_list(1:end-1)) , cos(angle_list(1:end-1)) ];

ArangeR = Amin(1):spacing:Amax(1);
ArangeG = Amin(2):spacing:Amax(2);
ArangeB = Amin(3):spacing:Amax(3);

Aall = generate_Avals(ArangeR, ArangeG);
[~, AvoteRG] = vote_2D(points(:,1:2), points_weight, directions_all, Aall, thres );

Aall = generate_Avals(ArangeG, ArangeB);
[~, AvoteGB] = vote_2D(points(:,2:3), points_weight, directions_all, Aall, thres );

Aall = generate_Avals(ArangeR, ArangeB);
[~, AvoteRB] = vote_2D(points(:,[1,3]), points_weight, directions_all, Aall, thres);

max_val = max( [max(AvoteRB(:)) , max(AvoteRG(:)) , max(AvoteGB(:)) ]);
AvoteRG2 = AvoteRG./max_val;
AvoteGB2 = AvoteGB./max_val;
AvoteRB2 = AvoteRB./max_val;
A11 = repmat( reshape(AvoteRG2, length(ArangeG),length(ArangeR))', 1,1,length(ArangeB));
tmp = reshape(AvoteRB2, length(ArangeB),length(ArangeR))';
A22 = repmat(reshape(tmp, length(ArangeR),1,length(ArangeB)) , 1,length(ArangeG),1);
tmp2 = reshape(AvoteGB2, length(ArangeB),length(ArangeG))';
A33 = repmat(reshape(tmp2, 1, length(ArangeG),length(ArangeB)) , length(ArangeR),1,1);
AvoteAll = A11.*A22.*A33;
[~, idx] = max(AvoteAll(:));
[idx_r,idx_g,idx_b] = ind2sub([length(ArangeR),length(ArangeG),length(ArangeB)],idx);
Aout = [ArangeR(idx_r), ArangeG(idx_g), ArangeB(idx_b)];


end 

function Aall = generate_Avals(Avals1, Avals2)

Avals1 = reshape(Avals1,[],1);
Avals2 = reshape(Avals2,[],1);
A1 = kron(Avals1, ones(length(Avals2),1));
A2 = kron(ones(length(Avals1),1), Avals2);
Aall = [A1, A2];
end 

function [Aout, Avote2] = vote_2D(points, points_weight, directions_all, Aall, thres)
n_directions = size(directions_all,1);
accumulator_votes_idx = false(size(Aall,1), size(points,1), n_directions);
for i_point = 1:size(points,1)
    for i_direction = 1:n_directions
        idx_to_use = find( (Aall(:, 1) > points(i_point, 1)) & (Aall(:, 2) > points(i_point, 2)));
        if isempty(idx_to_use), continue; end
		
        dist1 = sqrt(sum([Aall(idx_to_use, 1)-points(i_point, 1), Aall(idx_to_use, 2)-points(i_point, 2)].^2,2));
        dist1 = dist1./sqrt(2) + 1;
        
        dist =  -points(i_point, 1)*directions_all(i_direction,2) + ...
            points(i_point, 2)*directions_all(i_direction,1) + ...
            Aall(idx_to_use, 1)*directions_all(i_direction,2) - ...
            Aall(idx_to_use, 2)*directions_all(i_direction,1);
        idx = abs(dist)<2*thres.*dist1;
        if ~any(idx), continue; end

        idx_full = idx_to_use(idx);
        accumulator_votes_idx(idx_full, i_point,i_direction) = true;
    end
end
accumulator_votes_idx2 = (sum(uint8(accumulator_votes_idx),2))>=2; 
accumulator_votes_idx = bsxfun(@and, accumulator_votes_idx ,accumulator_votes_idx2);
accumulator_unique = zeros(size(Aall,1),1);
for iA = 1:size(Aall,1)
    idx_to_use = find(Aall(iA, 1) > points(:, 1) & (Aall(iA, 2) > points(:, 2)));
    points_dist = sqrt((Aall(iA,1) - points(idx_to_use,1)).^2+(Aall(iA,2) - points(idx_to_use,2)).^2);
    points_weight_dist = points_weight(idx_to_use).*(5.*exp(-reshape(points_dist,1,[]))+1); 
    accumulator_unique(iA) = sum(points_weight_dist(any(accumulator_votes_idx(iA,idx_to_use,:),3)));
end
[~, Aestimate_idx] = max(accumulator_unique);
Aout = Aall(Aestimate_idx,:);
Avote2 = accumulator_unique; 

end 