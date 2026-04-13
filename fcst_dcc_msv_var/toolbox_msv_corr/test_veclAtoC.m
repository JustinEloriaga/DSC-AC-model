
% Test code for the transformation

DD = randn(1,6);

tic 
for i=1:100
AA = veclAtoC_paper(DD);
end
toc

tic 
for i=1:100
BB = veclAtoC(DD);
end
toc

tic 
for i=1:100
CC = veclAtoC_4v_mex(DD);
end
toc