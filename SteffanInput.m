function C = SteffanInput(filename, headerlines, formatstring, filetype, delimiter);

%function to extract data from datafiles of various formats, one at a time
%FILENAME: the fullpath (or possibly relative path) of input file as a string.
%HEADERLINES: the number of lines that preceed actual data. These are either extracted as the header or entirely ignored.
%FORMATSTRING: used with textfiles to assign input/output styles to use for each column of data. e.g., formatstring='%n%n%n%n%n%s';
%FILETYPE: string which selects from filetype options
%C: output structure with C.data, C.words and C.header possible output variables. Exact format depends on which filetype used.

% TO DO
%1: want to have something that flexibly reads text files, rejects certain
%lines (e.g., header) and converts the rest of the information into
%sensible matrix (if all numbers), or possibly intellegently works out what
%is needed.
% 2.1: allow importing of eprime files
% 2.2: add functionality based on filename extension which allows importing of
% excel files (other types too?)
% 3: output format to be data.numbers and data.words ?

if exist('filetype','var')==0 %variable undefined?
    filetype='e-prime';
    %         filetype='text';
    
end
if exist('headerlines','var')==0 %variable undefined?
    headerlines=4;
    
end

switch filetype
    case 'excel'
        if exist('filename','var')==0 %variable undefined?
            filename='C:\Documents and Settings\skennett\My Documents\4Projects\Nick\TVcuts\Exp 2 xls files\P4.xls'; %a working file
        end
        %%% %%% %%% EXCEL %%% %%% %%%
        %This works for excel files, excluding headerlines (originally FROM SARTanalysis.m)
        %%MIGHT WORK FOR OTHER EXTENSIONS FORMATS THAT MATLAB CAN INTERNALLY COPE WITH...TRY IT AND SEE
        clear C
        h=importdata(filename);
        C.data=h.data(headerlines+1:end,:); %coded as numbers
        C.words=h.textdata; %not sure how this formatted
        C.header='Header not specifically extracted. See C.words';
        
    case 'text'
        if exist('filename','var')==0 %variable undefined?
            filename='C:\Documents and Settings\skennett\My Documents\4Projects\Rachel\Study1_raw.results\Discrimination_voicepilot1-2-1.edat'; %Rache;'s e-prime data
        end
        if exist('formatstring','var')==0 %variable undefined?
            formatstring='%s%n%n%n%n%s';
        end
        if exist('delimiter','var')==0 %variable undefined?
            delimiter=' ';
        end
        %%% %%% %%% TEXT OF ANY EXTENSION%%% %%% %%%
        clear C
        fid = fopen(filename);
        C.data=textscan(fid,formatstring,'Headerlines',headerlines, 'Delimiter', delimiter); %how to automate the input format string?
        C.header='Header ignored, not extracted';
        fclose(fid);
        
    case 'e-prime' %DOESN'T WORK IN EDAT FILES. only the e-prime text files
        if exist('filename','var')==0 %variable undefined?
            filename='C:\Documents and Settings\skennett\My Documents\4Projects\Rachel\Study1_raw.results\Discrimination_voicepilot1-2-1.txt'; %Rache;'s e-prime data
        end
        clear C
        C=importdata(filename);
end

debugvariable='finished'; %used to put a breakfpoint infront of when testing function and wanting to keep all variables in workspace.

%%% %%% %%% WORK IN PROGRESS %%% %%% %%%
% % % Take 2
% % % p= importdata(filename,'\t'); %all text from data file but not parsed!
% % [p pf]= importdata(filename); %all text from data file but not parsed!
% % % Take 3
% % [a,b,c,d] = finfo(filename);
% % q
%%%%FROM MARKUS BINDEMANN%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does this work for numberical tables/excel files?
% Does this work for any kind of text?
% is the "L=size(fline)" variable correct, shouldn't it be val not fline
% Why three calls within the loop?
% % % %fid = fopen(strcat(a,b,c));
% % % fid = fopen(filename);
% % % fline = fgets(fid);
% % % nbr = 1;
% % % while (fline ~= -1)
% % %     k = 1;
% % %     [val fline] = strtok(fline,char(9));
% % %     trialid(nbr,1) = cellstr(val);
% % %     [val fline] = strtok(fline,char(9));
% % %     trialid(nbr,2) = cellstr(val);
% % %     [val fline] = strtok(fline,char(9));
% % %     trialid(nbr,3) = cellstr(val);
% % %     L = size(fline);
% % %     while (L(2)>0)
% % %         [val fline] = strtok(fline,char(9));
% % %         emdata(nbr,k) = int32(str2num(val));
% % %         L = size(fline);
% % %         k = k+1;
% % %     end
% % %     fline = fgets(fid);
% % %     nbr = nbr + 1;
% % % end
% % % fclose(fid)