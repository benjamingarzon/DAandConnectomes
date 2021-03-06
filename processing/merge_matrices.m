function merge_matrices
  clear all
%
CONN_FOLDER='connectomes_noGSR' %'connectomes'
%subject_list='sub-D02 sub-D03 sub-D04 sub-D05 sub-D07 sub-D08 sub-D11 sub-D13 sub-D14 sub-D15 sub-D16 sub-D17 sub-D18 sub-D19 sub-D20 sub-D21 sub-D22 sub-D23 sub-D24 sub-D25 sub-D26 sub-D29 sub-D30 sub-D31 sub-D33 sub-D34 sub-D35 sub-D36 sub-D37 sub-D38 sub-D39 sub-D40 sub-D42 sub-D43 sub-D46 sub-D47 sub-D48 sub-D49 sub-D50 sub-D52 sub-D55 sub-D56 sub-D58 sub-D60 sub-D61 sub-D62 sub-D63 sub-D64 sub-D65 sub-D66 sub-D67 sub-D70 sub-D72 sub-D80 sub-D82 sub-D83 sub-D84 sub-D85 sub-D86 sub-D90';
matrix_file = 'zFC_150.csv';
merged_matrices_file = 'zFC_all_150.mat';
codes_file = '~/Data/DAD/parcellations/shen/Group_seg150_BAindexing_setA.txt';

subject_list = fullfile('~/Data/DAD/processed/fmriprep', CONN_FOLDER, 'Subjects_GNG.txt');
cd(fullfile('~/Data/DAD/processed/fmriprep/', CONN_FOLDER, 'GNG'));
merge_matrices_shen(subject_list, matrix_file, codes_file, merged_matrices_file);

subject_list = fullfile('~/Data/DAD/processed/fmriprep', CONN_FOLDER, 'Subjects_RS.txt');
cd(fullfile('~/Data/DAD/processed/fmriprep/', CONN_FOLDER, 'RS'));
merge_matrices_shen(subject_list, matrix_file, codes_file, merged_matrices_file);

subject_list = fullfile('~/Data/DAD/processed/fmriprep', CONN_FOLDER, 'Subjects_TAB.txt');
cd(fullfile('~/Data/DAD/processed/fmriprep/', CONN_FOLDER, 'TAB'));
merge_matrices_shen(subject_list, matrix_file, codes_file, merged_matrices_file);

end

function merge_matrices_shen(subject_list, matrix_file, codes_file, merged_matrices_file)
%function test(subject_list, matrix_file, codes_file, merged_matrices_file)

subjects = importdata(subject_list);

info = readtable(codes_file, 'Delimiter', ',');
labels = info.name;
regions = info.subunnit;

merged_matrices = [];
merged_matrices_mat = [];
for i=1:numel(subjects)

M = load(fullfile('.',subjects{i}, matrix_file));
M(isinf(M)) = 0;
if i==1
    merged_matrices = M;
else
    merged_matrices = cat(3, merged_matrices, M);
end
nans = sum(isnan(M), 2) == size(M, 1) - 1;
display([subjects{i}, ' ', num2str(sum(nans)) ] );

end
nans = any(isnan(merged_matrices), 3);
valid = sum(nans) < (numel(labels) - 1) ;
display(sum(valid))
valid_labels = labels(valid);
valid_regions = regions(valid);
valid_indices = find(valid);

for i=1:numel(subjects)    
    M = merged_matrices(:,:,i);
    M = M(valid, valid);
    M(isinf(M)) = 0;
    subjects{i} = subjects{i}(5:7);
    merged_matrices_mat = [merged_matrices_mat; squareform(M)];    
end
display(['Number of subjects: ' num2str(size(merged_matrices_mat, 1))])
display(['Total valid indices: ' num2str(numel(valid_labels))])
display(size(merged_matrices_mat))

save(merged_matrices_file, 'merged_matrices', 'labels', 'subjects','valid_regions','valid_indices', 'valid_labels', 'merged_matrices_mat') 

end