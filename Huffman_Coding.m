% This script uses and analyses Huffman codes to compress a text file.

% The first variant is the standard approach using a 0th order Markov 
% model, i.e. the code of any character being only dependent on its 
% frequency in the text. The second variant uses a 1st order Markov model 
% approach, i.e. the code of any character being only dependent on the 
% frequency after the character before.

% by Matthias Wolf - 2018


close all; clear; clc; tic;
fprintf('-------------------------------------------------------------\n');
fprintf('Huffman_Coding.m running...\n');
fprintf('-------------------------------------------------------------\n');


%% Prepare text file

% Read text file
text = fileread('testdata.txt');

% Split string into single characters
text_chars = split(text,{''});
text_chars = text_chars(2:end-1);
    % without this, whitespaces at the beginning and end are generated

% Define set of characters
char_set = {' ','0','1','2','3','4','5','6','7','8','9',...
    'A','B','C','D','E','F','G','H','I','J','K','L','M',...
    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

% Set of characters in uint8
uint8_set = double(cell2mat(char_set));
uint8_edges = [uint8_set-0.5,uint8_set(end)+0.5];
    % Bin edges: All values -0.5, last value +0.5

% Define category names
cat_set = {'blank','0','1','2','3','4','5','6','7','8','9',...
    'A','B','C','D','E','F','G','H','I','J','K','L','M',...
    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

% Generate categories
text_cat = categorical(text_chars,char_set,cat_set);


%% Huffman coding

% Compute entropy in case of i.i.d.
entr_iid = -log2(1/size(char_set,2));
fprintf('Entropy in case of i.i.d.: %.4f bits.\n',entr_iid);

% Compute histogram
hist = histcounts(text_cat);

% Plot histogram
figure('Name','Histogram',...
    'units','normalized','outerposition',[0 0 1 1]);
histogram('Categories',categorical(cat_set),'BinCounts',hist,...
    'FaceColor','k','FaceAlpha',1);
grid on;
xlabel('Character')
ylabel('Observations')
title('Histogram')

% Compute probabilities
prob = hist./sum(hist);

% Compute entropy
entropy = -sum(prob.*log2(prob));
fprintf('Entropy: %.4f bits.\n',entropy);

% Compute Huffman codes and average codeword length
[dict,avglen] = huffmandict(char_set,prob,2,'min');
fprintf('Average codeword length: %.4f bits.\n',avglen);

% Compute coding efficiency in percent
eff_pc = (entropy/avglen)*100;
fprintf('Efficiency: %.3f %%.\n',eff_pc);

% Encode text - CAUTION: Takes long time
encode = huffmanenco(text,dict);

% Decode text - CAUTION: Takes long time
decode = huffmandeco(encode,dict)';

% Determine if decoded text equals the original text
success = isequal(text_chars,decode);

% Print result
if success == 1 % in case of successful encoding and decoding
    fprintf('Successful encoding and decoding. The Huffman codes are:\n');
    for k = 1:size(char_set,2)
        fprintf('''%s'': ',cat_set{k});
        output=sprintf('%d',dict{k,2});
        fprintf('%s\n', output)
    end
else % in case of unseccessful encoding and decoding, throw error message
    fprintf('Error occured after encoding and decoding.\n')
end


%% 1st order Markov model

fprintf('-------------------------------------------------------------\n');
fprintf('1st order Markov model analysis:\n');

% Find pairs of letters that appear in the text
if mod(size(text_chars,1),2)==1 % odd numbers of characters
    two_chars_1 = reshape(text_chars(1:end-1),2,[])';
    two_chars_2 = reshape(text_chars(2:end),2,[])';
else % even number of characters
    two_chars_1 = reshape(text_chars(1:end),2,[])';
    two_chars_2 = reshape(text_chars(2:end-1),2,[])';
end
two_chars = [two_chars_1;two_chars_2];

% Since there is no histcounts2 support for 2D categorical data,
% convert all characters to its corresponding uint8 value
two_chars_uint8 = uint8(cell2mat(two_chars));

% Compute 2D histogram: Two following letters
hist2 = histcounts2(two_chars_uint8(:,1),two_chars_uint8(:,2),...
    uint8_edges,uint8_edges); % Bin edges

% Conditional probabilities:
% Probability of second character if first character appeared
cond_prob = hist2./sum(hist2,2);
test_sum = sum(cond_prob,2); % should be 1 everywhere if correct


%% Huffman codes for 1st order Markov model

% Pre-allocate Huffman codes and average codeword length for every case
markov_dict = cell(size(char_set,2),2);
markov_avglen = zeros(size(char_set,2),1);

% Compute Huffman codes and average codeword length for every case
for k = 1:size(char_set,2)
    markov_dict{k,1} = char_set{1,k}; % first column for first character
    [markov_dict{k,2},markov_avglen(k)] = ...
        huffmandict(char_set,cond_prob(k,:),2,'min');
        % second column for dictorinary in every case
end

% Compute entropy for every case
p_log_p = cond_prob.*log2(cond_prob);
p_log_p(~isfinite(p_log_p)) = 0;
markov_entropy = -sum(p_log_p,2);

% Difference between average codeword length and entropy in every case
markov_diff = markov_avglen - markov_entropy;

% Coding efficiency in percent in every case
markov_eff_pc = (markov_entropy./markov_avglen).*100;

% Compute averages of entropies, average codeword lengths and efficiencies
mean_avglen = mean(markov_avglen);
mean_entropy = mean(markov_entropy);
mean_eff_pc = (mean_entropy/mean_avglen)*100;

% Print average codeword length, entropy and efficiency in every case
fprintf('     | Average codeword length | '); % first line in table
fprintf('  Entropy   | Efficiency\n');        % first line in table
for k = 1:size(char_set,2) % table contents
    fprintf(' ''%s'' |       %.4f bits       | %.4f bits | %.3f %%\n',...
        char_set{k},markov_avglen(k),markov_entropy(k),markov_eff_pc(k));
end
fprintf('Mean |       %.4f bits       | %.4f bits | %.3f %%\n',...
    mean_avglen,mean_entropy,mean_eff_pc); % prints the mean values


%% Plots for 1st order Markov model

% Plot 2D histogram - CAUTION: Takes long time
figure('Name','1st order Markov model: 2D Histogram',...
    'units','normalized','outerposition',[0 0 1 1]); 
set(gcf,'renderer','zbuffer'); % for better performance
scatterbar3((size(hist2,1):-1:1).*ones(size(hist2)),...
    (size(hist2,2):-1:1)'.*ones(size(hist2)),hist2,0.8);
    % external function, see file exchange doc
colormap('jet');
colorbar;
axis tight;
grid on;
xticks(1:size(cat_set,2));
xticklabels(fliplr(cat_set)); % flipped for better view
xlabel('Second character');
yticks(1:size(cat_set,2));
yticklabels(fliplr(cat_set)); % flipped for better view
ylabel('First character');
zlabel('Observations');
title('1st order Markov model: 2D Histogram');

% Plot conditional probabilities: 3D bar plot - CAUTION: Takes long time
figure('Name','1st order Markov model: Conditional probabilities',...
    'units','normalized','outerposition',[0 0 1 1]); 
set(gcf,'renderer','zbuffer'); % for better performance
scatterbar3((size(cond_prob,1):-1:1).*ones(size(cond_prob)),...
    (size(cond_prob,2):-1:1)'.*ones(size(cond_prob)),cond_prob,0.8);
    % external function, see file exchange doc
colormap('jet');
colorbar;
axis tight;
grid on;
xticks(1:size(cat_set,2));
xticklabels(fliplr(cat_set)); % flipped for better view
xlabel('Second character');
yticks(1:size(cat_set,2));
yticklabels(fliplr(cat_set)); % flipped for better view
ylabel('First character');
zlabel('Conditional probability');
title('1st order Markov model: Conditional probabilities');

% Plot conditional probabilities: stacked bars - CAUTION: Takes long time
figure('Name','1st order Markov model: Conditional probabilities',...
    'units','normalized','outerposition',[0 0 1 1]); 
bar(cond_prob,'stacked'); % stacked bars should add up to value 1
axis tight;
colormap('jet');
legend(cat_set,'Location','westoutside');
xticks(1:size(cat_set,2));
xticklabels(cat_set);
xlabel('First character');
ylabel('Conditional probability of second character');
title('1st order Markov model: Conditional probabilities');


%% Execution time

execution_time = toc;
fprintf('-------------------------------------------------------------\n');
fprintf('Execution time: %.0f seconds.\n',execution_time);
fprintf('-------------------------------------------------------------\n');
fprintf('\n\n');

% END OF SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
