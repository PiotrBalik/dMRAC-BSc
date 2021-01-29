function teks = arr2lat(X,name)
[row,l]=size(X);
% teks=[n '=\blkmat{' sprintf(['%2.4G ' repmat('& %2.4G ',1,l-1) '\\\\ \n'],X') '}\\'];
styl={'mrac','pid'};
P1=['\begin{table}[h!]\centering\caption{}\begin{tabular}{|c|' repmat('c|',1,l+1) '}'];
teks=[sprintf(['& & a%d & %.2G ' repmat('& %.2G ',1,l-2) '\\\\ \\hline \n'],[[mod(2:row+1,2)+1]'  X(:,2:end)]') ];
teks=strrep(teks,'a1','mrac');
teks=strrep(teks,'a2','pid');
P2=['\end{tabular}\label{tab:' name '}\end{table}'];

ff=fopen([pwd '\wyniki\tabele.txt'],'a');

teks=[P1 teks P2];

fwrite(ff,teks);
fclose(ff);
% teks=sprintf(['%2.4G ' repmat('& %2.4G ',1,l-1) '\\\\ \n'],X');
% clipboard('copy',teks);
end