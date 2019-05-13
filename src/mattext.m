function [text]=mattex(A);

[n,m]=size(A);

text='\\begin{table}[htbp] \n';
texttmp='\\centering';
text=strcat(text,texttmp);
texttmp='\\begin{tabular}{|';
text=strcat(text,texttmp);
    for j=1:m+1
        texttmp='c|';
        text=strcat(text,texttmp);
    end
texttmp='}\n';
text=strcat(text,texttmp);
for i=1:n
    texttmp=' \\hline \t';
    text=strcat(text,texttmp);
        texttmp='$';
        text=strcat(text,texttmp);
    texttmp=num2str(i);
    text=strcat(text,texttmp);
            texttmp='$';
        text=strcat(text,texttmp);
    texttmp='  \t&';
    text=strcat(text,texttmp);
    for j=1:m
        if j>1 
            texttmp=' \t & \\;';
            text=strcat(text,texttmp);
        end    
        texttmp='$';
        text=strcat(text,texttmp);
        texttmp=num2str(A(i,j));
        text=strcat(text,texttmp);
        texttmp='$';
        text=strcat(text,texttmp);

    end
    texttmp=' \\\\ \n';
    text=strcat(text,texttmp);
end
texttmp='\\hline \n \\end{tabular} \n';
text=strcat(text,texttmp);
texttmp='\\end{table} ';
text=strcat(text,texttmp);

sprintf(text)



