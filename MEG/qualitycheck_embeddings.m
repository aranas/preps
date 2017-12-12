% This is a script to visualize and test relationships in the MOUS word
% embeddings

load /project/3011020.09/MEG/misc/mous_stimuli.mat
indsel      = [];
vec         = [];
word        = {};
for i = 1:length(tlck_sent.trialinfo)
    wavid   = tlck_sent.trialinfo(i,6);
    posword = tlck_sent.trialinfo(i,5);
    label   = strtok(stimuli(wavid).words(posword).POS,'(');
    vec     = stimuli(wavid).words(posword).word2vec;
    if ~isempty(vec) && strcmp(label{1},'N')
        word2v(i,:)     = vec;
        indsel          = [indsel i];
        word{i}         = stimuli(wavid).words(posword).word;
    end
end
word2v    = word2v(indsel,:);
word      = word(indsel);


% results in 236 words
% do dimensionality reduction
% needed tsne function for this
% I used this toolbox https://lvdmaaten.github.io/drtoolbox/code/drtoolbox.tar.gz
addpath(genpath('/home/language/sopara/matlab/drtoolbox'))

word2v_r     = compute_mapping(word2v,'tSNE',3);
figure;
plot3(word2v_r(:,1),word2v_r(:,2),word2v_r(:,3),'.')
text(word2v_r(:,1),word2v_r(:,2),word2v_r(:,3),[word{:}])

%cluster words using kmeans to see if groupings are meaningful
[IDX, C] = kmeans(word2v, 20);
groups = unique(IDX);
colors = hsv(length(groups));
% Initialize some axes


% Plot each group individually:
for k = 1:length(groups)
    view(3)
grid on
hold on
    % Get indices of this particular unique group:
    ind = IDX==groups(k);
    % Plot only this group:
    plot3(word2v_r(ind,1),word2v_r(ind,2),word2v_r(ind,3),'.','color',colors(k,:),'markersize',20);
    text(word2v_r(ind,1),word2v_r(ind,2),word2v_r(ind,3),[word{ind}],'FontSize',20);
    pause
end


