function f = build_cloth_features(map,sbin)
%BUILD_CLOTH_FEATURES
if isempty(map)
    f = [];
    return;
end

% List of clothing id
cloth_ids = [0 5 7 9 12 15 18 22 24 26 27 31 36 49 53 58 63 64 74 88 ...
    90 99 101 115 117 154 159 173 174 190 193 229 241 245 251 296 308 ...
    318 323 336 341 526 643 678 723 775 814 969 973 1082 1086 1145 ...
    1226 1577 22628 22629];
%cloth_ids = internal_cloth_id;
%map = internal_cloth_id(map);

% Split map into blocks
nc = repmat(sbin,1,floor((size(map,2)-1)/sbin+1));
nr = repmat(sbin,1,floor((size(map,1)-1)/sbin+1));
nc(end) = mod(size(map,2)-1,sbin)+1;
nr(end) = mod(size(map,1)-1,sbin)+1;
map_c = mat2cell(map,nr,nc);

% Calculate histograms for each
B = zeros(numel(map_c),numel(cloth_ids));
for i = 1:numel(map_c)
    B(i,:) = histc(map_c{i}(:),cloth_ids);
end
B = reshape(B,[size(map_c),numel(cloth_ids)]);

% Sum up neighboring 4 blocks
out_h = max(round(size(map,1)/sbin)-2,0);
out_w = max(round(size(map,2)/sbin)-2,0);
f = zeros(out_h,out_w,numel(cloth_ids));
for j = 1:size(f,1)
    for i = 1:size(f,2)
        h = B(j,i,:)+B(j+1,i,:)+B(j,i+1,:)+B(j+1,i+1,:);
        %h = cat(3,B(j,i,:),B(j+1,i,:),B(j,i+1,:),B(j+1,i+1,:));
        sh = sum(h(:));
        assert(sh>0);
        f(j,i,:) = h(:)./sh;
    end
end

end

%%
function Y = internal_cloth_id(X)
% convert fashionista clothing id to internal id
% see below for the corresponding internal id

persistent M;
if isempty(M)
    % garment id
    cloth_ids = ...
        [117 15 26 1145 173 241 90 643 27 74 88 58 53 251 678 190 ...
        159 973 526 323 1577 1086 723 18 969 336 22628 318 814 341 1082 775 ...
        99 31 1226 49 193 64 392 12 101 174 229 5 923 7 24 0 9 245 115 22 ...
        1401 36 308 63 296 154 22629]+1; % +1 since it contains 0
    % group them together
    group_ids = [1 2 2 2 3 4 4 4 5 6 6 6 7 7 7 8 8 9 9 9 9 9 9 9 9 9 10 11 ...
        11 11 11 12 12 13 13 13 13 13 13 13 13 14 14 14 15 15 15 16 17 17 ...
        17 17 17 17 17 17 17 17 18];
    M = sparse(cloth_ids,1,group_ids);
end

if nargin < 1
    Y = full(unique(M));
else
    Y = full(M(X+1));
end

end

% Internal-name Fashionista-name Fashionista-id Internal-id
% accessories	accessories	117	1
% bag	bag	15	2
% bag	purse	26	2
% bag	wallet	1145	2
% belt	belt	173	3
% body-accessory	necklace	241	4
% body-accessory	scarf	90	4
% body-accessory	tie	643	4
% boots	boots	27	5
% bottomware	jeans	74	6
% bottomware	leggings	88	6
% bottomware	pants	58	6
% dress	dress	53	7
% dress	jumper	251	7
% dress	romper	678	7
% eyeware	glasses	190	8
% eyeware	sunglasses	159	8
% footware	clogs	973	9
% footware	flats	526	9
% footware	heels	323	9
% footware	loafers	1577	9
% footware	pumps	1086	9
% footware	sandals	723	9
% footware	shoes	18	9
% footware	sneakers	969	9
% footware	wedges	336	9
% hair	hair	22628	10
% hand-accessory	bracelet	318	11
% hand-accessory	gloves	814	11
% hand-accessory	ring	341	11
% hand-accessory	watch	1082	11
% head-acceessory	earrings	775	12
% head-acceessory	hat	99	12
% inner	blouse	31	13
% inner	bodysuit	1226	13
% inner	bra	49	13
% inner	intimate	193	13
% inner	shirt	64	13
% inner	swimwear	392	13
% inner	t-shirt	12	13
% inner	top	101	13
% legware	socks	174	14
% legware	stockings	229	14
% legware	tights	5	14
% lower-body	panties	923	15
% lower-body	shorts	7	15
% lower-body	skirt	24	15
% null	null	0	16
% outer	blazer	9	17
% outer	cape	245	17
% outer	cardigan	115	17
% outer	coat	22	17
% outer	hoodie	1401	17
% outer	jacket	36	17
% outer	suit	308	17
% outer	sweater	63	17
% outer	sweatshirt	296	17
% outer	vest	154	17
% skin	skin	22629	18
