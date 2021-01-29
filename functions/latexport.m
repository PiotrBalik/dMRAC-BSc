function latexport(hfull,titles)

% wrap1='\begin{figure*}\centering%\resizebox{0.90\linewidth}{!}{}}\includegraphics[width=\textwidth]{';
wrap1='\begin{figure}[h!]\centering\includegraphics[width=\textwidth]{';
wrap2='}\label{fig:';
wrap3='}\caption{';
wrap4='}\end{figure}';
len=length(hfull);
texts=cell(1,2*len);
for figs=1:len
    texts{2*figs-1}=' ';
    texts{2*figs} = [wrap1 sprintf([ 'Rysunki/' titles{figs} '%d.png'],figs)...
          wrap2 sprintf([titles{figs} '%d'],figs) wrap3 titles{figs} wrap4];

    text = [pwd sprintf([ '\\wykresy\\' titles{figs} '%d.png'],figs)];
    print(text,hfull(figs),'-dpng')
end
ff=fopen([pwd '\wyniki\T' titles{1} '.txt'],'a');
fwrite(ff,cell2mat(texts));
fclose(ff);
% clipboard('copy',cell2mat(texts))
end