function [text]=mattex(A,texto);

[n,m]=size(A);

text='\\[ \n';
texttmp='=\\left[\n';
text=strcat(text,texttmp);
texttmp='\\begin{array}{';
text=strcat(text,texttmp);
    for j=1:m
        texttmp='c';
        text=strcat(text,texttmp);
    end
texttmp='}\n';
text=strcat(text,texttmp);
for i=1:n
    for j=1:m
        if j>1
            texttmp=' \t& ';
            text=strcat(text,texttmp);
        end    
        texttmp=num2str(A(i,j));
        text=strcat(text,texttmp);

    end
    texttmp=' \\\\ \n';
    text=strcat(text,texttmp);
end
texttmp='\\end{array} \n';
text=strcat(text,texttmp);
texttmp='\\right] \n';
text=strcat(text,texttmp);
texttmp='\\] ';
text=strcat(text,texttmp);

sprintf(text)

