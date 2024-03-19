function [ outputStruct ] = AutoKEGGRec( organismCodes, varargin )
% 
% Assembles one or more genome-scale metabolic reconstructions COBRA
% structures based on KEGG which are output to workspace and can be
% written to file. Reconstructions can be for single organisms,
% communities, or consolidated, i.e. the union of several organisms.
% 
% USAGE:
% 
%   outputStruct = AutoKEGGRec(organismKEGGIDs, varargin)
% 
% INPUTS:
%
%   organismKEGGIDs: string array of KEGG organism IDs
% 
% OPTIONAL INPUTS:
% 
%   varargin:   Optional flags, altering the output structure.
%               - 'ConsolidatedRec': Creates a consolidated
%                 reconstruction of the given organisms.
%               - 'SingleRecs': Creates separate reconstructions for all
%                 given organisms.
%               - 'CommunityRec': Creates a community reconstruction for
%                 all given organisms.
%               - 'writeSBML': Uses the COBRA writeSBML function to save the
%                 reconstructions as SBML files in your "Current Folder".
%               - 'OmittedData': Adds data to the output structure on
%                 reactions omitted during reconstruction.
%               - 'OrgRxnGen': Adds the Organisms-Reactions-Genes-Matrix for the
%                 query organisms to the output structure.
%               - 'DisconnectedReactions': Adds a cell to the output
%                 structure containing the KEGG reaction IDs for all
%                 reactions in the metabolic network not connected to the
%                 giant component.
%               - 'GenePlot': Plots the number of reactions associated with
%                 X number of genes against X.
%               - 'Histogram': Plots number of reactions present in X
%                 organisms against X.
% 
% .. AUTHORS:
%   - Emil Karlsen and Christian Schulz, April 2018 - Created
%   - Christian Schulz, August 2018 - Added compound assembly and annotations + notes + changes in omitted output
%   - Emil karlsen, March 2024 - Compatibility patch for KEGG API updates + various minor changes
% 
% EXAMPLES:
% 
% % Create a consolidated reconstruction of the E. coli K12 strains in KEGG
% % using most of the flags
% output = AutoKEGGRec(["eco","ecj","ecd","ebw","ecok"], 'ConsolidatedRec', 'SingleRecs', 'CommunityRec', 'OrgRxnGen', 'GenePlot', 'Histogram', 'OmittedData', 'DisconnectedReactions')
% 
% % Create both a consolidated and a community reconstruction of the E.
% % coli K12 strains in KEGG
% output =
% AutoKEGGRec(["eco","ecj","ecd","ebw","ecok"],'ConsolidatedRec','CommunityRec')
% 
% % Create separate reconstructions of the E. coli K12 strains in KEGG
% output = AutoKEGGRec(["eco","ecj","ecd","ebw","ecok"],'SingleRecs')
% 
% NOTE:
% 
%   Organism IDs have to be KEGG IDs (e.g. "eco" or "T00007").
%   They have to be given in a string array (e.g. ["eco","ecj","ecd"]).
%
%   Refer to the Manual for further and more detailed instructions and to the paper for general explanations.
% 
%   New versions may be available:
% 
%           https://www.ntnu.edu/almaaslab/downloads
%           https://github.com/emikar/AutoKEGGRec
%
%   
%
%%

disp('This version of AutoKEGGRec is the published version, and might not be the most recent.');
disp('Please check the websites specified within the help function for the newest version.');
pause(5)
tic

%% Input args

ConsolidatedRecFlag=false;
SingleRecsFlag=false;
CommunityRecFlag=false;
writeSBMLflag=false;

OrgRxnGenFlag = false;
GenePlotFlag = false;
HistogramFlag = false;
OmittedDataFlag = false;
DisconnectedReactionsFlag = false;
DebugFlag = false;

arguments=false;
recBeingBuilt = false;
% Check OrganismCodes input
if ~isstring(organismCodes)
   error('The given Organism Codes have to be a list of strings!')
end
if length(char(organismCodes))<3
   error('Missing or wrong organism codes.')
end

KeyWords = {'ConsolidatedRec', 'SingleRecs', 'CommunityRec', 'writeSBML', 'OrgRxnGen', 'GenePlot', 'Histogram', 'OmittedData', 'DisconnectedReactions', 'Debug'};
%Check flags

if numel(varargin) > 0
   for argum=1:numel(varargin)
       if ischar(varargin{argum})
           if any(ismember(varargin{argum}, KeyWords))
               arguments = true;            
           else
               disp(varargin{argum})
               error('Unrecognized command: "%s"\nPlease refer to the help function (F1) or the manual.',string(varargin{argum}));
           end
       end
   end
else
   arguments = true;
end



if arguments
    if isempty(varargin)
        disp("Default settings: Consolidated reconstruction as output.")
        ConsolidatedRecFlag = true;
        recBeingBuilt = true;
    else
        disp("The selected flags:")
        for argum=1:numel(varargin)
            n = find(strcmp(KeyWords,varargin(argum)));
            switch n
                case 1
                    ConsolidatedRecFlag = true;
                    recBeingBuilt = true;
                    disp("ConsolidatedRecflag")
                case 2
                    SingleRecsFlag = true;
                    recBeingBuilt = true;
                    disp("SingleRecsflag")
                case 3
                    CommunityRecFlag = true;
                    recBeingBuilt = true;
                    disp("CommunityRecflag")
                case 4
                    writeSBMLflag = true;
                    disp("writeSBMLflag")
                case 5
                    OrgRxnGenFlag = true;
                    disp("OrgRxnGenFlag")
                case 6
                    GenePlotFlag = true;
                    disp("GenePlotFlag")
                case 7
                    HistogramFlag = true;
                    disp("HistogramFlag")
                case 8
                    OmittedDataFlag = true;
                    disp("OmittedDataFlag")
                case 9
                    DisconnectedReactionsFlag = true;
                    disp("DisconnectedReactionsFlag")
                case 10
                    DebugFlag = true;
                    disp("DebugFlag")
                    % Causes the annotation matrices to be written to file for inspection.
                    % Particularly useful to check for changes in the KEGG API.
            end
        end
    end
else
    error('Bad arguments');
end

if (writeSBMLflag) && ~recBeingBuilt
    error('In order to write an SBML file, a reconstruction needs to be built.\nUse optional input ConsolidatedRec, SingleRecs, and/or CommunityRec. %s',"")
elseif (DisconnectedReactionsFlag) && ~recBeingBuilt
    error('In order to identify disconnected reactions, a reconstruction needs to be built.\nUse optional input ConsolidatedRec, SingleRecs, and/or CommunityRec. %s',"")
end

%% Initialization (defining variables etc.)

apiRequestWaitTime = 0.1;

options = weboptions('Timeout', 500);
keggURL = "http://rest.kegg.jp/";

%Building URLs for retrieving organism gene and enzyme data
nOrganisms = length(organismCodes);
organismGenesToECURLs = strings(1,nOrganisms);

fprintf("\n \n");
disp("KEGG organism codes:")
disp(organismCodes)

for i=1:nOrganisms
    organismGenesToECURLs(i) = keggURL+"link/"+organismCodes(i)+"/ec";
end
fprintf("\n \n");
disp("Starting accessing the KEGG Database for the requested organisms.")

%Building URLs for retrieving enzyme-reaction and reaction-metabolite
%linkage data.
ecToRxnURL = keggURL+"link/ec/rn";
rxnToMetURL = keggURL+"link/rn/cpd";

% Retrieving reaction, enzyme, and compound lists
rxnList = cellstr(webread(keggURL+"list/reaction", options));
enzList = cellstr(webread(keggURL+"list/enzyme", options));
metList = cellstr(webread(keggURL+"list/compound", options));

% Splitting into lines 
rxnList = splitlines(rxnList);
enzList = splitlines(enzList);
metList = splitlines(metList);

% Splitting lines by tab
rxnList = split(rxnList, '	');
enzList = split(enzList, '	');
metList = split(metList, '	');

% Keeping first element of each line
rxnList = rxnList(:,1);
enzList = enzList(:,1);
metList = metList(:,1);

% Adding prefixes to the KEGG IDs
rxnList = strcat("rn:", rxnList);
enzList = strcat("ec:", enzList);
metList = strcat("cpd:", metList);

nRxns = length(rxnList);
nEnzs = length(enzList);
nMets = length(metList);


%Building reaction, enzyme and compound references (name to #)
rxnRef = containers.Map(rxnList,1:length(rxnList));
enzRef = containers.Map(enzList,1:length(enzList));
metRef = containers.Map(metList,1:length(metList));

disp("All organism URLs have been built, and the reactions, compound and EC number linkage have been downloaded.")

%% Makes a cell array of string matrices associating EC numbers and genes

organismGeneToECLists = cell(1,nOrganisms);
ecToRxn = cell(1);
rxnToMet = cell(1);

% Make an nRxns-by-2 matrix of strings with EC numbers and associated genes
% for each organism, and retrieve enzyme-reaction and reaction-metabolite linkage

for i=1:nOrganisms+2
    if i<=nOrganisms
        urlIn = organismGenesToECURLs(i);
    elseif i==nOrganisms+1
        urlIn = ecToRxnURL;
    elseif i==nOrganisms+2
        urlIn = rxnToMetURL;
    end
    importedRawData = webread(urlIn, options);
    splitData = strsplit(importedRawData);
    spliDatLen = length(splitData);
    if mod(spliDatLen,2)==1
        iterNum = spliDatLen-1;
    else
        iterNum = spliDatLen;
    end
    linkageMat = repmat("",[iterNum/2 2]);
    for j=1:iterNum
        linkageMat(ceil(j/2),mod(j+1,2)+1) = splitData(j);
    end
    tableThing = table(linkageMat(:,1),linkageMat(:,2));
    if i<=nOrganisms
        organismGeneToECLists{i} = table2cell(tableThing);
    elseif i==nOrganisms+1
        ecToRxn = table2cell(tableThing);
    elseif i==nOrganisms+2
        rxnToMet = table2cell(tableThing);
    end
end
time = toc;
fprintf("\nTime so far running AutoKEGGRec: %.0f seconds.\n\n", time);
disp("All organism gene-to-reaction linkages have been downloaded.")

%Building reaction-enyme matrix
rxnEnzMat = zeros(length(rxnList),length(enzList));
for i=1:length(ecToRxn)
    rxnPos = rxnRef(char(ecToRxn{i,1}));
    enzPos = enzRef(char(ecToRxn{i,2}));
    rxnEnzMat(rxnPos,enzPos)=1;
end

%Building reaction-compound matrix
rxnMetMat = zeros(length(rxnList),length(metList));
for i=1:length(rxnToMet)
    metPos = metRef(char(rxnToMet{i,1}));
    rxnPos = rxnRef(char(rxnToMet{i,2}));
    rxnMetMat(rxnPos,metPos)=1;
end

disp("Compound-to-reaction and reaction-to-enzyme matrices built. Starting organism-to-reaction matrix now.")

%% Building reaction-organism matrix
% Builds a binary matrix tying each reaction in KEGG to each of the query
% organisms

rxnOrganismMatrix = zeros(nOrganisms,nRxns);

for i=1:nOrganisms
    for j=1:length(organismGeneToECLists{i})
        reactions = find(rxnEnzMat(:,enzRef(char(cellstr(organismGeneToECLists{i}(j,1)))))>0);
        for k=1:length(reactions)
            rxnOrganismMatrix(i,rxnRef(rxnList{reactions(k)})) = 1;
        end
    end
end

disp("Matrices built.")

if recBeingBuilt
    disp("Commencing assembly of reconstruction(s).");
end

% Selecting all reactions in at least one organism
if nOrganisms == 1
   rxnsToBuildList = find(rxnOrganismMatrix>0);
else
   rxnsToBuildList = find(sum(rxnOrganismMatrix([1:end],:)>0));
end

%% Building reaction-organism matrix with genes

% Builds a string matrix tying each reaction in KEGG to each of the query
% organisms through the respective genes.
% 
% NB: The reactions are assumed to depend on the genes in an OR-fashion
% NB2: ATP reaction
% 
% Ex:
%         Org1        Org2        ...
% Rx1     G1 | G2     G14 | G10
% Rx2     -           G12
% :

rxnOrganismGeneMatrix = strings(nRxns,nOrganisms);
rxnOrganismGeneMatrixCounter = zeros(nRxns,nOrganisms);

for org=1:nOrganisms
    for gene=1:length(organismGeneToECLists{org})
        ecNum = char(string(organismGeneToECLists{org}(gene,1)));
        ecIndex = enzRef(ecNum);
        rxnIndicesForEC = find(rxnEnzMat(:,ecIndex)>0);
        for rxnIndex=1:length(rxnIndicesForEC)
            if rxnOrganismGeneMatrix(rxnIndicesForEC(rxnIndex),org)~=""
                rxnOrganismGeneMatrix(rxnIndicesForEC(rxnIndex),org) = rxnOrganismGeneMatrix(rxnIndicesForEC(rxnIndex),org)+" | ";
            end
            geneName = string(organismGeneToECLists{org}(gene,2));
            geneName = strrep(geneName,organismCodes(org)+":","");
            %genes without organism code HERE
            rxnOrganismGeneMatrix(rxnIndicesForEC(rxnIndex),org) = rxnOrganismGeneMatrix(rxnIndicesForEC(rxnIndex),org)+geneName;
            rxnOrganismGeneMatrixCounter(rxnIndicesForEC(rxnIndex),org) = rxnOrganismGeneMatrixCounter(rxnIndicesForEC(rxnIndex),org)+1;
        end
    end
end

%% Plot number of genes per reaction for each organism

if GenePlotFlag
    geneCounterList = zeros(max(max(rxnOrganismGeneMatrixCounter)),nOrganisms);
    reactionCounter = 0;
    subplotX = ceil(sqrt(nOrganisms));
    subplotY = ceil(nOrganisms/ceil(sqrt(nOrganisms)));
    figure
    for org=1:nOrganisms
        for rxn=1:nRxns
           if rxnOrganismGeneMatrixCounter(rxn,org)>0
               geneCounterList(rxnOrganismGeneMatrixCounter(rxn,org),org)=geneCounterList(rxnOrganismGeneMatrixCounter(rxn,org),org)+1;
               reactionCounter = reactionCounter+1;
           end
        end
        subplot(subplotX,subplotY,org)
        geneCountsToPlot = geneCounterList;
        geneCountsToPlot(max(find(geneCountsToPlot(:,org)~=0)),org)=0;
        geneCountToPlot = geneCountsToPlot(1:max(find(geneCountsToPlot(:,org)~=0)),org);
        bar(geneCountToPlot)
        xlabel(['Number of genes per reaction']) % x-axis label
        ylabel('Number of reactions') % y-axis label
        set(gca, 'YScale', 'log')
        title(['KEGG organism ID: ' organismCodes(org)])
    end
end

%% Plot rxns for organism histogram

if HistogramFlag
    figure;
    organismHistData = sum(rxnOrganismMatrix);
    histogram(organismHistData(organismHistData>0));
    title(string('Number of Organisms (KEGG IDs: ' + strjoin(organismCodes,', ') + ') sharing the same reaction'))
    xlabel("Number of organisms sharing number of reactions")
    ylabel("Number of reactions")
    histogramPlotHandle = gca;
    set(histogramPlotHandle, 'YScale', 'log')
    xLims = xlim;
    set(histogramPlotHandle, 'xtick', 1:ceil(xLims(2)))
end

%% Getting all the reactions from KEGG according to matrix

if recBeingBuilt || OmittedDataFlag
    
    fprintf("%-55s %10d\n", "Number of KEGG reactions for the KEGG organism(s) queried:", length(rxnsToBuildList));
    disp("Downloading all reactions from KEGG (this takes a while).")

    fprintf('Progress:\n');
    fprintf(['\n' repmat('.',1,100) '\n\n']);
    loadBarUpdatePoints = ceil(linspace(1,length(rxnsToBuildList)));

    downloadedData = cell(length(rxnsToBuildList));
    for i=1:length(rxnsToBuildList)
        downloadedData{i} = cellstr(webread(keggURL+"get/"+rxnList(rxnsToBuildList(i)), options));
        downloadedData{i} = replace(downloadedData{i},"///","");
        pause(apiRequestWaitTime)
        if ~isempty(find(loadBarUpdatePoints==i, 1))
            fprintf('\b|\n');
        end
    end

    equation_list_all = strings(1,length(rxnsToBuildList));
    metaDataList = equation_list_all;
    metaDataDumpList = equation_list_all;
    reactionDumpList = equation_list_all;

    glycanReactions = 0;
    nContainingReactions = 0;
    numGenReactions = 0;
    numMovedReactions = 0;
    sugarReactions = 0;
    polymerReactions = 0;

    for i=1:length(downloadedData)
        webRxnData = string(downloadedData{i});
        webRxnData = splitlines(webRxnData);
        generalReac = webRxnData(contains(webRxnData,"COMMENT"));
        rxnEquation = webRxnData(contains(webRxnData,"EQUATION"));
        if ~isempty(webRxnData(contains(webRxnData,"REMARK")))
            numMovedReactions = numMovedReactions + 1;
            if ~isempty(rxnEquation(contains(rxnEquation, 'G')))
                reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
                metaDataDumpList(i) = join(string(webRxnData),";");
                metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction contains glycans, which have a C-compound ID as well. Check carefully. Check carefully!"]);
                glycanReactions = glycanReactions+1;
                % Removing reactions containing sugars (double reactions)
            elseif ~isempty(rxnEquation(contains(rxnEquation, 'n')))
                reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
                metaDataDumpList(i) = join(string(webRxnData),";");
                metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction is a reaction containing n in the equation. Check carefully!"]);
                nContainingReactions = nContainingReactions+1;
                % Removing reactions containing polymers (double reactions)
            else
                equation_list_all(i) = replace(rxnEquation,"EQUATION    ","");
                numMovedReactions = numMovedReactions - 1;
                metaDataList(i) = join(string(webRxnData),";");
                %This way all the "same as suggar" , ... are sorted out, only
                %the C-compounts are kept.
            end
        elseif ~isempty(generalReac(contains(generalReac,'GENERIC')))
            numGenReactions = numGenReactions + 1;
            reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
            metaDataDumpList(i) = join(string(webRxnData),";");
            metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction is a generic reaction. Check carefully!"]);
            %Generic reactions are removed
        elseif ~isempty(generalReac(contains(generalReac,'GENERAL')))
            numGenReactions = numGenReactions + 1;
            reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
            metaDataDumpList(i) = join(string(webRxnData),";");
            metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction is a genral reaction. Check carefully!"]);
            %Generic reactions are removed
        elseif ~isempty(generalReac(contains(generalReac,'general reaction')))
            numGenReactions = numGenReactions + 1;
            reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
            metaDataDumpList(i) = join(string(webRxnData),";");
            metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction is a genral reaction. Check carefully!"]);
            %Generic reactions are removed
        else
            if ~isempty(rxnEquation(contains(rxnEquation, 'G')))
                % Reactions containing suggars that are not marked within the
                % comment field; They are added to the reconstruction, but they are
                % marked within the attention field.
                sugarReactions = sugarReactions + 1;
                equation_list_all(i) = replace(rxnEquation,"EQUATION    ","");
                metaDataList(i) = join(string(webRxnData),";");
                metaDataList(i) = insertBefore(metaDataList(i),[";DEFINITION"],[";ATTENTION    The reaction contains a sugar!"]);
            elseif ~isempty(rxnEquation(contains(rxnEquation, 'n')))
                % Reactions containing "n" within the reaction equation. This
                % is meant to deal with polymer reactions. These are removed,
                % since for the S-Matrix A + B <-> C + A becomes a problem. The
                % user has the possibility to implement them later.
                %
                % First, the reactions are changed:
                % nA + (n-1)B <-> C + (n+1)A becomes
                % A + B <-> C + A
                % The original equation is stored in the comments.
                polymerReactions = polymerReactions + 1;
                rxnEquation1 = string(rxnEquation);
                rxnEquation1 = strsplit(rxnEquation1);
                EQ = ("EQUATION    ");
                for ii=1:length(rxnEquation1)
                    totest = char(rxnEquation1(ii));
                    if ~isempty(totest(contains(totest, 'C'))) || ~isempty(totest(contains(totest, 'G')))
                        if ~isempty(totest(contains(totest, 'C')))
                            posC = strfind(totest, "C");
                            if ii<length(rxnEquation1)
                                trytest = rxnEquation1(ii+1);
                                if ~isempty(trytest(contains(trytest, '=')))
                                    EQ = EQ + totest(posC:posC+5);
                                else
                                    temp = totest(posC:posC+5) + " + ";
                                    EQ = EQ + temp;
                                end
                            else
                            EQ = EQ + totest(posC:posC+5);  
                            end

                        elseif ~isempty(totest(contains(totest, 'G')))
                            posG = strfind(totest, "G");
                            if ii<length(rxnEquation1)
                                trytest = rxnEquation1(ii+1);
                                if ~isempty(trytest(contains(trytest, '=')))
                                    EQ = EQ + totest(posG:posG+5);
                                else
                                    temp = totest(posG:posG+5) + " + ";
                                    EQ = EQ + temp;
                                end
                            else
                            EQ = EQ + totest(posG:posG+5);  
                            end
                        end
                    elseif ~isempty(totest(contains(totest, '=')))
                        totest = " " + totest + " ";
                        EQ = EQ + totest;    
                    end
                end
                rxnEquation = EQ;
                reactionDumpList(i) = replace(rxnEquation,"EQUATION    ","");
                metaDataDumpList(i) = join(string(webRxnData),";");
                metaDataDumpList(i) = insertBefore(metaDataDumpList(i),[";DEFINITION"],[";ATTENTION    The reaction is a polymer reaction. For the reconstruction, the polymerisation has been removed. The original reaction is stored!"]);
            else
                % "Normal" reactions are put into the List
                equation_list_all(i) = replace(rxnEquation,"EQUATION    ","");
                rxnEquation = equation_list_all(i);
                metaDataList(i) = join(string(webRxnData),";");
                %Checking for equatic aberrations, ie. characters in
                %biochemical equations that were not accounted for above.
                charsToTest = char(rxnEquation);
                greenLightedChars = '1234567890C<=>+ ';
                for k=1:length(charsToTest)
                    if ~contains(greenLightedChars,charsToTest(k))
                        errorMsg = sprintf("Error using AutoKEGGRec: Unexpected sign in reaction equaton: %s\n",rxnEquation);
                        error(char(errorMsg));
                    end
                end
            end
        end
        if ~isempty(find(loadBarUpdatePoints==i))
            fprintf('\b|\n');
        end
    end
    
    fprintf("\n \n")
    time = toc;
    disp("All reactions needed have been fetched.")
    fprintf("\nTime so far running AutoKEGGRec: %.0f seconds.\n\n", time);
    
    %Reaction checking for double Compounds (e.g. C00029 + C00760 <=> C00015 +
    %C00760); Moving them into the reactionDumpList
    reactionsProducingTheirOwnSubstrate = 0;
    for rxn=1:length(equation_list_all)
        onlycompounds = false;
        cpdNames = split(equation_list_all(rxn));   
        for i=1:length(cpdNames)
            if length(char(cpdNames(i)))==6 && sum(count(cpdNames, cpdNames(i)))>1
                        onlycompounds = true;
            end
        end
        if onlycompounds
            reactionDumpList(rxn) = equation_list_all(rxn);
            equation_list_all(rxn) = "";
            metaDataDumpList(rxn) = metaDataList(rxn);
            metaDataList(rxn) = "";
            reactionsProducingTheirOwnSubstrate = reactionsProducingTheirOwnSubstrate + 1;
            metaDataDumpList(rxn) = insertBefore(metaDataDumpList(rxn),[";DEFINITION"],[";ATTENTION    The reaction hase the same compound on both sites of the reaction. Check carefully!"]);
        end
    end

    rxnKeggNamesList = string(rxnList(rxnsToBuildList));


    %% Removing the empty fields, that only the reactions stay within the matrix to not check for compounds again

    deletionIndex = zeros(1,length(equation_list_all));
    equation_list_keep = equation_list_all;
    for i=1:length(equation_list_all)
        if length(char(equation_list_keep(i)))<1
            deletionIndex(i) = 1;
        end
    end

    % Removing reactions not to be kept
    equation_list_keep(find(deletionIndex)) = [];

    %% Getting the Compound information from KEGG
    
    % Getting the compounds
    cpdlist = [];
    for reac=1:length(equation_list_keep)
        reactemp = split(equation_list_keep(reac));
        for cmpd=1:length(reactemp)
            if length(char(reactemp(cmpd))) == 6
                cpdlist = [cpdlist, reactemp(cmpd)];
            end
        end
    end
    cpdlist = unique(cpdlist);
    
    fprintf("%-55s %10d\n", "Number of unque KEGG compounds for the organism(s) queried:", length(cpdlist));
    disp("Downloading all compounds from KEGG (this takes a while).")
    
    % downloading the compound information from KEGG
    downloadedDataCPD = cell(size(cpdlist));
    fprintf('Progress: \n');
    fprintf(['\n' repmat('.',1,100) '\n\n']);
    loadBarUpdatePoints = ceil(linspace(1,length(cpdlist)));
    for i=1:length(cpdlist)
        downloadedDataCPD{i} = cellstr(webread(keggURL+"/get/"+cpdlist(i), options));
        downloadedDataCPD{i} = replace(downloadedDataCPD{i},"///","");
        pause(apiRequestWaitTime)
        if ~isempty(find(loadBarUpdatePoints==i))
            fprintf('\b|\n');
        end
    end

    stringlisttosearchforCPD = {'ENTRY', 'NAME', 'FORMULA', 'MASS', 'MOL_WEIGHT', 'REMARK','COMMENT', 'REACTION', 'PATHWAY', 'MODULE', 'ENZYME', 'BRITE', 'DBLINKS', 'ATOM', 'BOND', 'COMPOSITION', 'NODE', 'EDGE'};
    
    fprintf("\n \n")
    time = toc;
    disp("All compounds needed have been fetched.")
    fprintf("\nTime so far running AutoKEGGRec: %.0f seconds.\n\n", time);
    
    disp("Adding KEGG annotations to the compounds.")
    CPDannotationMatrix = num2cell(zeros(length(downloadedDataCPD),length(stringlisttosearchforCPD)));

    % Preparing the Annotation fields by reshaping the KEGG information into
    % a cell of specified fields according to stringlisttosearchforCPD
    for i=1:length(downloadedDataCPD)
        WebCPDData = string(downloadedDataCPD{i});
        WebCPDData = join(split(WebCPDData));
        containingmodulesCPD=[];
        for list=1:length(stringlisttosearchforCPD)
            searchvarCPD = strfind(WebCPDData, stringlisttosearchforCPD(list));
            if ~isempty(searchvarCPD)
                containingmodulesCPD = [containingmodulesCPD, string(stringlisttosearchforCPD(list))];
            end
        end
        annotationsStrListCPD = WebCPDData;
        for inputelementinstringCPD=1:length(containingmodulesCPD)
            annotationsStrListCPD = insertBefore(annotationsStrListCPD, containingmodulesCPD(inputelementinstringCPD), ["%$%"]);
            annotationsStrListCPD = insertAfter(annotationsStrListCPD, containingmodulesCPD(inputelementinstringCPD), ["%$%"]);
        end
        annotationsStrListCPD=strsplit(annotationsStrListCPD, "%$%");
        for inputelementinstringCPD=1:length(annotationsStrListCPD)
            position_of_annotation = strcmp(stringlisttosearchforCPD, annotationsStrListCPD(inputelementinstringCPD));
            if sum(position_of_annotation) == 1
                CPDannotationMatrix(i, position_of_annotation) = {annotationsStrListCPD(inputelementinstringCPD+1)};
            end
        end
    end

    % Changing the compound name fields to the corresponding KEGG IDs
    for cpd=1:length(CPDannotationMatrix(:,1))
        splitstringCPD = split(string(CPDannotationMatrix(cpd)));
        CPDannotationMatrix(cpd) = cellstr(splitstringCPD(2));
    end
    
    % Identify generic compounds according to their molelcular weight = 0
    genericCPDs = string(CPDannotationMatrix(find(string(CPDannotationMatrix(1:end,4))=="0")));

    % Identify reactions containing these generic compounds
    reactionCheckList = [];
    for reac=1:length(equation_list_keep)
        reactemp = split(equation_list_keep(reac));
        for cmpd=1:length(reactemp)
            if length(char(reactemp(cmpd))) == 6 && contains(reactemp(cmpd), genericCPDs)
                reactionCheckList = [reactionCheckList, equation_list_keep(reac)];
            end
        end
    end
    reactionCheckList = unique (reactionCheckList);
    % Cleaning the reaction list of these reactions and adding them to the
    % omitted reactions, adding a comment to the cleaned reactions, why
    % they have been moved
    for cleanreac = 1:length(equation_list_all)
        if sum(strcmp(equation_list_all(cleanreac), reactionCheckList)) == 1
           reactionDumpList(cleanreac) =  string(equation_list_all(cleanreac));
           equation_list_all(cleanreac) = "";
           metaDataDumpList(cleanreac) = metaDataList(cleanreac);
           if contains(metaDataDumpList(cleanreac), "ATTENTION")
               metaDataDumpList(cleanreac) = insertBefore(metaDataDumpList(cleanreac),[";DEFINITION"],[" Furthermore, this reaction contains a generic compound or a compound without mass in KEGG. Check the compounds carefully before adding reaction to model!"]);
           else
               metaDataDumpList(cleanreac) = insertBefore(metaDataDumpList(cleanreac),[";DEFINITION"],[";ATTENTION    This reactions contains a generic compound or a compound without mass in KEGG. Check the compounds carefully before adding reaction to model!"]);
           end
           metaDataList(cleanreac) = "";
        end
    end

    % Removing reactions not to be kept from the reaction list
    deletionIndex = zeros(1,length(equation_list_all));
    equation_list_keep = equation_list_all;
    for i=1:length(equation_list_all)
        if length(char(equation_list_keep(i)))<1
            deletionIndex(i) = 1;
        end
    end
    equation_list_keep(find(deletionIndex)) = [];

    %% Removing the empty fields, that only the reactions stay within the matrix
    deletionIndexDUMP = zeros(1,length(reactionDumpList));
    reactionDumpList2 = reactionDumpList;
    for i=1:length(reactionDumpList)
        if length(char(reactionDumpList2(i)))<1
            deletionIndexDUMP(i) = 1;
        end
    end
    reactionDumpList2(find(deletionIndexDUMP)) = [];

    if OmittedDataFlag
        dumpAnnotations=cellstr(metaDataDumpList);
        dumpAnnotations = dumpAnnotations(~cellfun('isempty',dumpAnnotations));
    end
end

%% Adding annotations to the omitted output
if OmittedDataFlag
    %%% Compounds
    omittedCPDs = struct();
    for dumpcpds=1:length(genericCPDs)
        omittedCPDs.(char(genericCPDs(dumpcpds)))= cell(length(stringlisttosearchforCPD),2);
        posinCPDList = strcmp(CPDannotationMatrix(:,1), genericCPDs(dumpcpds));
        for downw=1:length(stringlisttosearchforCPD)
            omittedCPDs.(char(genericCPDs(dumpcpds))){downw,1} = string(stringlisttosearchforCPD(downw));
            omittedCPDs.(char(genericCPDs(dumpcpds))){downw,2} = string(CPDannotationMatrix(posinCPDList,downw));
        end
    end
    
    %%% Reactions
    stringlisttosearchfor = {'ENTRY', 'NAME', 'ATTENTION', 'DEFINITION', 'EQUATION', 'COMMENT', 'RCLASS', 'ENZYME', 'PATHWAY', 'ORTHOLOGY', 'DBLINKS', 'REFERENCE', 'RHEA:'};
    disp("Adding KEGG annotations to the reactions and compounds within the omitted data.")
    for line=1:length(dumpAnnotations)
        annotationsStrListdump=string(dumpAnnotations(line));
        containingmodules={};
        for list=1:length(stringlisttosearchfor)
            searchvar = strfind(annotationsStrListdump, stringlisttosearchfor(list));
            if ~isempty(searchvar)
                containingmodules{end+1} = string(stringlisttosearchfor(list));
            end
        end
        containingmodules=string(containingmodules);
        annotationssplitsting=annotationsStrListdump;
        for inputelementinstring=1:length(containingmodules)
            if inputelementinstring==1
                annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(2) + "  "));
                relevantannotation = annotationssplitsting(1);
                relevantannotation=strsplit(relevantannotation, (containingmodules(inputelementinstring)));
                relevantannotation=relevantannotation(2);
                temp2=relevantannotation;
                temp1='';
                while ~strcmp(temp1,temp2)
                    temp1=temp2;
                    temp2=regexprep(temp1,'  ',' ');
                end
                temp2=char(strtrim(temp2));
                relevantannotation=temp2(1:6);
                annotationssplitsting = annotationssplitsting(2);
            elseif inputelementinstring==length(containingmodules)
                if containingmodules(inputelementinstring)=='RHEA:'
                    annotationssplitsting=strsplit(annotationssplitsting, (":"));
                    relevantannotation = annotationssplitsting(2);
                    temp2=relevantannotation;
                    temp1='';
                    while ~strcmp(temp1,temp2)
                        temp1=temp2;
                        temp2=regexprep(temp1,'  ',' ');
                    end
                    temp2=strtrim(temp2);
                    relevantannotation=temp2;
                else
                    annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(inputelementinstring) + " "));
                    relevantannotation = annotationssplitsting(1);
                    temp2=relevantannotation;
                    temp1='';
                    while ~strcmp(temp1,temp2)
                        temp1=temp2;
                        temp2=regexprep(temp1,'  ',' ');
                    end
                    temp2=strtrim(temp2);
                    relevantannotation=temp2;
                end
            else
                tempsplitnumber = inputelementinstring + 1;
                annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(tempsplitnumber) + "  "));
                relevantannotation = annotationssplitsting(1);
                temp2=relevantannotation;
                temp1='';
                while ~strcmp(temp1,temp2)
                    temp1=temp2;
                    temp2=regexprep(temp1,'  ',' ');
                end
                temp2=strtrim(temp2);
                relevantannotation=temp2;
                if length(annotationssplitsting)>1
                    annotationssplitsting = annotationssplitsting(2);
                else
                    annotationssplitsting=relevantannotation;
                    relevantannotation = 0;
                end
            end
            position_of_annotation = strcmp(stringlisttosearchfor, containingmodules(inputelementinstring));
            annotationMatrixdump(line, position_of_annotation) = {relevantannotation};
        end 
    end
    omittedrxns = struct();
    for dumprxns=1:length(annotationMatrixdump(:,2))
        omittedrxns.(char(annotationMatrixdump(dumprxns,1))) = cell(length(stringlisttosearchfor),2);
        for downw=1:length(stringlisttosearchfor)
            omittedrxns.(char(annotationMatrixdump(dumprxns,1))){downw,1} = string(stringlisttosearchfor(downw));
            if isempty(annotationMatrixdump{dumprxns, downw})
                omittedrxns.(char(annotationMatrixdump(dumprxns,1))){downw,2}= [];
            else
                omittedrxns.(char(annotationMatrixdump(dumprxns,1))){downw,2} = string(annotationMatrixdump(dumprxns, downw));
            end
        end
    end
end

pause(3)

%% Building consolidated reconstruction
if recBeingBuilt
    disp("Building consolidated reconstruction.");
    pause(2)
    annotations=cellstr(reshape(metaDataList,[length(metaDataList),1]));

    % Removing empty annotations
    annotations = annotations(~cellfun('isempty',annotations));
    
    ReactionFormulas = cellstr(equation_list_keep);
    ReactionNames = cellstr(rxnKeggNamesList);
    ReactionNames2 = cellfun(@(x) x(1:end,4:9), ReactionNames, 'un',0);
    ReactionNames2(find(deletionIndex)) = [];
    lowerbounds = [ones(1,length(equation_list_keep))*-20];
    upperbounds = [ones(1,length(equation_list_keep))*20];
    [throwAwayText,consolidatedStruct] = evalc(char("createModel(ReactionNames2, ReactionNames2, ReactionFormulas, 'lowerBoundList', lowerbounds, 'upperBoundList', upperbounds);"));

    %% Adding annotations to reactions
    stringlisttosearchfor = {'ENTRY', 'NAME', 'ATTENTION', 'DEFINITION', 'EQUATION', 'COMMENT', 'RCLASS', 'ENZYME', 'PATHWAY', 'MODULE', 'BRITE', 'ORTHOLOGY', 'DBLINKS', 'REFERENCE', 'RHEA:'};
    annotationsStrList=strings(1,length(annotations));
    annotationMatrix = num2cell(zeros(length(annotations),length(stringlisttosearchfor)));

    disp("Adding KEGG annotations to the reactions.")
    for line=1:length(annotations)
        annotationsStrList=string(annotations(line));
        containingmodules={};
        for list=1:length(stringlisttosearchfor)
            searchvar = strfind(annotationsStrList, stringlisttosearchfor(list));
            if ~isempty(searchvar)
                containingmodules{end+1} = string(stringlisttosearchfor(list));
            end
        end
        containingmodules=string(containingmodules);
        annotationssplitsting=annotationsStrList;
        for inputelementinstring=1:length(containingmodules)
            if inputelementinstring==1
                annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(2) + "  "));
                relevantannotation = annotationssplitsting(1);
                relevantannotation=strsplit(relevantannotation, (containingmodules(inputelementinstring)));
                relevantannotation=relevantannotation(2);
                temp2=relevantannotation;
                temp1='';
                while ~strcmp(temp1,temp2)
                    temp1=temp2;
                    temp2=regexprep(temp1,'  ',' ');
                end
                temp2=char(strtrim(temp2));
                relevantannotation=temp2(1:6);
                annotationssplitsting = annotationssplitsting(2);
            elseif inputelementinstring==length(containingmodules)
                if containingmodules(inputelementinstring)=='RHEA:'
                    annotationssplitsting=strsplit(annotationssplitsting, (":"));
                    relevantannotation = annotationssplitsting(2);
                    temp2=relevantannotation;
                    temp1='';
                    while ~strcmp(temp1,temp2)
                        temp1=temp2;
                        temp2=regexprep(temp1,'  ',' ');
                    end
                    temp2=strtrim(temp2);
                    relevantannotation=temp2;
                else
                    annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(inputelementinstring) + " "));
                    relevantannotation = annotationssplitsting(1);
                    temp2=relevantannotation;
                    temp1='';
                    while ~strcmp(temp1,temp2)
                        temp1=temp2;
                        temp2=regexprep(temp1,'  ',' ');
                    end
                    temp2=strtrim(temp2);
                    relevantannotation=temp2;
                end
            else
                tempsplitnumber = inputelementinstring + 1;
                annotationssplitsting=strsplit(annotationssplitsting, (";" + containingmodules(tempsplitnumber) + "  "));
                relevantannotation = annotationssplitsting(1);
                temp2=relevantannotation;
                temp1='';
                while ~strcmp(temp1,temp2)
                    temp1=temp2;
                    temp2=regexprep(temp1,'  ',' ');
                end
                temp2=strtrim(temp2);
                relevantannotation=temp2;
                if length(annotationssplitsting)>1
                    annotationssplitsting = annotationssplitsting(2);
                else
                    annotationssplitsting=relevantannotation;
                    relevantannotation = 0;
                end
            end
            position_of_annotation = strcmp(stringlisttosearchfor, containingmodules(inputelementinstring));
            annotationMatrix(line, position_of_annotation) = {relevantannotation};
        end 
    end

    for rxn=1:length(annotationMatrix(:,2))
        if string(annotationMatrix(rxn,2)) == "0"
            annotationMatrix(rxn,2) = cellstr("Name not available. KEGG ID: "+string(cellstr(string(annotationMatrix(rxn,1)))));
        end
    end

    % Writing annotationMatrix to file for debug
    if DebugFlag
       writetable(cell2table(annotationMatrix), 'annotationMatrix.csv', 'Delimiter', '\t')
    end
    
    % COBRA supported fields according to https://github.com/opencobra/cobratoolbox/blob/master/docs/source/notes/COBRAModelFields.md
    dynamic_index = strcmp(stringlisttosearchfor, "ENTRY");
	consolidatedStruct.rxnKeggID = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "NAME");
	consolidatedStruct.rxnNames = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "REFERENCE");
    consolidatedStruct.rxnReferences = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "ENZYME");
    consolidatedStruct.rxnECNumbers = cellstr(string(annotationMatrix(:,dynamic_index)));

    dynamic_index = strcmp(stringlisttosearchfor, "PATHWAY");
    consolidatedStruct.subSystems = cell(length(cellstr(string(annotationMatrix(:,dynamic_index)))),1);
    for rxn=1:length(cellstr(string(annotationMatrix(:,dynamic_index))))
        consolidatedStruct.subSystems{rxn} = cellstr(string(annotationMatrix(rxn,dynamic_index)));
    end

    %COBRA not supported fields
    dynamic_index = strcmp(stringlisttosearchfor, "ATTENTION");
    consolidatedStruct.rxnAttention = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "DEFINITION");
    consolidatedStruct.rxnDefinition = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "EQUATION");
    consolidatedStruct.rxnEquation = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "COMMENT");
    consolidatedStruct.rxnComment = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "RCLASS");
    consolidatedStruct.rxnRclass = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "MODULE");
    consolidatedStruct.rxnRhea = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "BRITE");
    consolidatedStruct.rxnRhea = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "ORTHOLOGY");
    consolidatedStruct.rxnOrthology = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "DBLINKS");
    consolidatedStruct.rxnDBLinks = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "RHEA:");
    consolidatedStruct.rxnRhea = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "BRITE");
    consolidatedStruct.rxnBRITE = cellstr(string(annotationMatrix(:,dynamic_index)));
    dynamic_index = strcmp(stringlisttosearchfor, "MODULE");
    consolidatedStruct.rxnModule = cellstr(string(annotationMatrix(:,dynamic_index)));

    %% Adding annotations to the KEGG compounds

    % COBRA supported fields according to https://github.com/opencobra/cobratoolbox/blob/master/docs/source/notes/COBRAModelFields.md
    consolidatedStruct.metKeggID = cellfun(@(x) x(1:end-3), consolidatedStruct.mets, 'un', 0);
    consolidatedStruct.metFormulas = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metNames = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metNotes = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metChEBIID = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metPubChemID = cell(size(consolidatedStruct.mets,1),1);

    consolidatedStruct.metFormulas(:) = cellstr("");
    consolidatedStruct.metNames(:) = cellstr("");
    consolidatedStruct.metNotes(:) = cellstr("");
    consolidatedStruct.metChEBIID(:) = cellstr("");
    consolidatedStruct.metPubChemID(:) = cellstr("");

    % COBRA not supported fields
    consolidatedStruct.metMass = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metMolWeight = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metComment = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metReactionsContainingCPD = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metPathwaysContainingCPD = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metModulesContainingCPD = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metEnzymesContainingCPD = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metBRITE = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metDBLINKS = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metComposition = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metStructure2D.atomPosition = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metStructure2D.Bonds = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metStructure2D.node = cell(size(consolidatedStruct.mets,1),1);
    consolidatedStruct.metStructure2D.edge = cell(size(consolidatedStruct.mets,1),1);
    
    if DebugFlag
        % write CPDannotationMatrix to file
        writetable(cell2table(CPDannotationMatrix), 'CPDannotationMatrix.csv', 'Delimiter', '\t')
    end

    for cpds=1:length(consolidatedStruct.mets)
        if sum(contains(string(CPDannotationMatrix(:, 1)),string(consolidatedStruct.metKeggID(cpds))))==1
            positionofCPD = strcmp(string(CPDannotationMatrix(:, 1)),string(consolidatedStruct.metKeggID(cpds)));
            for cpdannoline=1:length(CPDannotationMatrix(positionofCPD, 1:end))
                switch cpdannoline
                    case 1
                        %KEGG ID
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metKeggID(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metKeggID(cpds) = cellstr("Unknown");
                        end
                    case 2
                        %Names
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metNames(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metNames(cpds) = cellstr("Unknown");
                        end
                    case 3
                        %Formula
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metFormulas(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metFormulas(cpds) = cellstr("Unknown");
                        end
                    case 4
                        %MASS or EXACT_MASS, depending on Compound or
                        %Glycan
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metMass(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            % This should never happen
                            fprintf("This message should never appear while adding compound annotation to the model! Compound: %s, Massfield: %s (Should not be zero)\n", string(consolidatedStruct.metKeggID(cpds)), string(CPDannotationMatrix(positionofCPD,cpdannoline)));
                            fprintf("Within the compound %s, something went wrong, since it appears to have no mass but is included into the model.... CHECK! \n", string(CPDannotationMatrix(positionofCPD,cpdannoline)));
                            consolidatedStruct.metExactMass(cpds) = cellstr("Unknown - Unexpected behavior");
                        end
                    case 5
                        %MOL_WEIGHT
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metMolWeight(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metMolWeight(cpds) = cellstr("Unknown");
                        end
                    case 6
                        %REMARK becomes Notes in COBRA
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metNotes(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metNotes(cpds) = cellstr("Unknown");
                        end
                    case 7
                        %COMMENT
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metComment(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metComment(cpds) = cellstr("Unknown");
                        end
                    case 8
                        %REACTION
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metReactionsContainingCPD(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metReactionsContainingCPD(cpds) = cellstr("Unknown");
                        end
                    case 9
                        %PATHWAY
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metPathwaysContainingCPD(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metPathwaysContainingCPD(cpds) = cellstr("Unknown");
                        end
                    case 10
                        %MODULE
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metModulesContainingCPD(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metModulesContainingCPD(cpds) = cellstr("Unknown");
                        end
                    case 11
                        %ENZYME
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metEnzymesContainingCPD(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metEnzymesContainingCPD(cpds) = cellstr("Unknown");
                        end
                    case 12
                        %BRITE
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metBRITE(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metBRITE(cpds) = cellstr("Unknown");
                        end
                    case 13
                        %DBLINKS
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metDBLINKS(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                            if contains(string(CPDannotationMatrix(positionofCPD,cpdannoline)), "PubChem") || contains(string(CPDannotationMatrix(positionofCPD,cpdannoline)), "ChEBI")
                                linesplitDBLINKS = split(string(CPDannotationMatrix(positionofCPD,cpdannoline)));
                                if contains(string(CPDannotationMatrix(positionofCPD,cpdannoline)), "PubChem")
                                    posPubChem = find(strcmp(linesplitDBLINKS, "PubChem:"));
                                    consolidatedStruct.metPubChemID(cpds) = cellstr(linesplitDBLINKS(posPubChem+1));
                                end
                                if contains(string(CPDannotationMatrix(positionofCPD,cpdannoline)), "ChEBI")
                                    posChEBIID = find(strcmp(linesplitDBLINKS, "ChEBI:"));
                                    consolidatedStruct.metChEBIID(cpds) = cellstr(linesplitDBLINKS(posChEBIID+1));
                                end
                            end
                        else
                            consolidatedStruct.metDBLINKS(cpds) = cellstr("Unknown");
                        end
                    case 14
                        %ATOM
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metStructure2D.atomPosition(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metStructure2D.atomPosition(cpds) = cellstr("Unknown");
                        end
                    case 15
                        %BOND
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metStructure2D.Bonds(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metStructure2D.Bonds(cpds) = cellstr("Unknown");
                        end
                    case 16
                        % COMPOSITION
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metComposition(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metComposition(cpds) = cellstr("Unknown");
                        end
                    case 17
                        %NODE
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metStructure2D.node(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metStructure2D.node(cpds) = cellstr("Unknown");
                        end
                    case 18
                        %EDGE
                        if length(char(string(CPDannotationMatrix(positionofCPD,cpdannoline)))) > 1
                            consolidatedStruct.metStructure2D.edge(cpds) = cellstr(CPDannotationMatrix(positionofCPD,cpdannoline));
                        else
                            consolidatedStruct.metStructure2D.edge(cpds) = cellstr("Unknown");
                        end
                    otherwise
                        error('New field within the KEGG compound annotations. Update the fields in the corresponding loop within the pipeline');
                end
            end
        else
            fprintf("This KEGG ID %s seem not to exist within the annotations downloaded from KEGG. Something went wrong, please check \n", string(consolidatedStruct.metKeggID(cpds)));
        end
    end

    %% Finding which reactions remain for each organism

    singleOrganismRxns = cell(nOrganisms,1);
    for org=1:nOrganisms
        relRxns = false(nRxns,1);
        relRxns(contains(rxnList,consolidatedStruct.rxns)) = true;
        relRxns(find((sum(rxnOrganismMatrix)>0)-(rxnOrganismMatrix(org,:)>0)>0)) = false;
        relInds = find(relRxns>0);
        singleOrganismRxns{org} = relInds;
    end

    %% Gene and protein annotations
    % ie. "model.genes", "model.grRules", "model.geneNames",
    % "model.proteinNames" and "model.proteins"

    % Declaring different struct subfields for gene and protein annotations

    consolidatedStruct.geneAnnotations.genes = cell(nOrganisms,1);
    consolidatedStruct.geneAnnotations.grRules = cell(nOrganisms,1);
    consolidatedStruct.geneAnnotations.geneNames = cell(nOrganisms,1);

    consolidatedStruct.proteinAnnotations.proteinNames = cell(nOrganisms,1);
    consolidatedStruct.proteinAnnotations.proteins = cell(nOrganisms,1);

    % Adding gene annotations and grRules

    for org=1:nOrganisms
        rxnCounter = 0;
        for rxn = 1:length(singleOrganismRxns{org})
            if rxnOrganismGeneMatrix(singleOrganismRxns{org}(rxn),org)~=""
                rxnCounter = rxnCounter+1;
            end
        end
    
        % Adding grRules
        consolidatedStruct.geneAnnotations.grRules{org} = strings(rxnCounter,1);
        for rxn=1:rxnCounter
            consolidatedStruct.geneAnnotations.grRules{org}(rxn) = rxnOrganismGeneMatrix(singleOrganismRxns{org}(rxn),org);
        end
    
        % Adding genes; not all genes are necessarily part of grRules
        consolidatedStruct.geneAnnotations.genes{org} = strings(length(unique(string(organismGeneToECLists{org}(:,2)))),1);
        consolidatedStruct.geneAnnotations.geneNames{org} = strings(length(unique(string(organismGeneToECLists{org}(:,2)))),1);

        consolidatedStruct.geneAnnotations.genes{org} = unique(string(organismGeneToECLists{org}(:,2)));
        consolidatedStruct.geneAnnotations.geneNames{org} = unique(string(organismGeneToECLists{org}(:,2)));
    end

    % Downloading gene-protein linkage

    organismGeneProteinData = cell(nOrganisms,1);

    for org=1:nOrganisms
        geneProtLink = strsplit(webread("http://rest.kegg.jp/list/"+organismCodes(org), options),'\n');
        if string(geneProtLink(end))==""
            geneProtLink = geneProtLink(1:length(geneProtLink)-1);
        end
        nGenesTot = size(geneProtLink,2);
        geneProteinList = strings(nGenesTot,2);
        for gene=1:nGenesTot
            geneProtEntries = strsplit(string(geneProtLink(gene)),'\t');
            geneProteinList(gene,1) = geneProtEntries(1);
            geneProteinList(gene,2) = geneProtEntries(2);
        end
        organismGeneProteinData{org} = geneProteinList;
    end

    % Adding protein annotation and updating gene names
    fprintf("Adding protein annotations and updating gene names for each organism: \n");
    for org=1:nOrganisms        
        consolidatedStruct.proteinAnnotations.proteinNames{org} = strings(length(unique(string(organismGeneToECLists{org}(:,2)))),1);
        consolidatedStruct.proteinAnnotations.proteins{org} = strings(length(unique(string(organismGeneToECLists{org}(:,2)))),1);
        fprintf("\t Organism: %s \n", string(organismCodes(org)));
        for gene=1:length(unique(string(organismGeneToECLists{org}(:,2))))
            geneAndProteinNames = organismGeneProteinData{org}(find(organismGeneProteinData{org}(:,1)==consolidatedStruct.geneAnnotations.genes{org}(gene)),2);
            if contains(geneAndProteinNames,";")
                geneAndProteinNames = strsplit(geneAndProteinNames,"; ");
                consolidatedStruct.geneAnnotations.geneNames{org}(gene) = geneAndProteinNames(1);
                consolidatedStruct.proteinAnnotations.proteinNames{org}(gene) = geneAndProteinNames(2);
                consolidatedStruct.proteinAnnotations.proteins{org}(gene) = geneAndProteinNames(2);
            else
                consolidatedStruct.proteinAnnotations.proteinNames{org}(gene) = geneAndProteinNames;
                consolidatedStruct.proteinAnnotations.proteins{org}(gene) = geneAndProteinNames;
            end
        end
    end
end
pause(3)

%% Building single reconstructions

if SingleRecsFlag || CommunityRecFlag
    disp("Building single reconstructions")
    singleStructs = struct();
    for org=1:nOrganisms
        orgName = organismCodes(org);
        % Removing reactions
        singleOrganismStruct = struct(consolidatedStruct);
        rxnIndicesForRemoval = false(nRxns,1);
        rxnIndicesForRemoval(contains(rxnList,consolidatedStruct.rxns)) = true;
        rxnIndicesForRemoval(singleOrganismRxns{org}) = false;
        rxnIndicesForRemoval = find(rxnIndicesForRemoval>0);
        rxnsToRemove = strings(length(rxnIndicesForRemoval),1);
        for rxn=1:length(rxnIndicesForRemoval)
            rxnToRemove = rxnList(rxnIndicesForRemoval(rxn));
            rxnToRemove = replace(string(rxnToRemove),"rn:","");
            rxnsToRemove(rxn) = rxnToRemove;
        end
        temp1 = zeros(3,1);
        temp1(1) = length(singleOrganismStruct.rxns);
        rxnsToRemove = cellstr(rxnsToRemove);
        temp1(2) = length(rxnsToRemove);
        singleOrganismStruct = removeRxns(singleOrganismStruct,rxnsToRemove);
        temp1(3) = length(singleOrganismStruct.rxns);

        % Adding gene and protein annotations
        singleOrganismStruct.genes = cellstr(consolidatedStruct.geneAnnotations.genes{org});
        singleOrganismStruct.genes = cellfun(@(x) replace(x,orgName+":",""), singleOrganismStruct.genes, 'un',0); % Removing organism code from gene names
        singleOrganismStruct.grRules = cellstr(consolidatedStruct.geneAnnotations.grRules{org});
        singleOrganismStruct.geneNames = cellstr(consolidatedStruct.geneAnnotations.geneNames{org});
        singleOrganismStruct.proteinNames = cellstr(consolidatedStruct.proteinAnnotations.proteinNames{org});
        singleOrganismStruct.proteins = cellstr(consolidatedStruct.proteinAnnotations.proteins{org});
        % Removing redundant fields
        singleOrganismStruct = rmfield(singleOrganismStruct,'geneAnnotations');
        singleOrganismStruct = rmfield(singleOrganismStruct,'proteinAnnotations');
        singleOrganismStruct = rmfield(singleOrganismStruct,'rxnGeneMat');
        % Adding structs to collection
        singleStructs.(char(orgName)) = singleOrganismStruct;
    end
    fieldNames = fieldnames(singleStructs);
end

%% Building community reconstruction

if CommunityRecFlag
    disp("Building community reconstruction");
    pause(2);
    %Building community reconstruction
    communityRxnNum = 0;
    communityMetNum = 0;
    communityGenNum = 0;
    for org=1:nOrganisms
        communityRxnNum = communityRxnNum+length(singleStructs.(fieldNames{org}).rxns);
        communityMetNum = communityMetNum+length(singleStructs.(fieldNames{org}).mets);
        communityGenNum = communityGenNum+length(singleStructs.(fieldNames{org}).genes);
    end
    communityRxnNames = cell(1,communityRxnNum);
    communityRxnEquations = cell(communityRxnNum,1);
    counter = 1;
    for org=1:nOrganisms
        codeToAdd = organismCodes(org);
        updateInterval = counter:(counter-1+length(singleOrganismRxns{org}));
        rxnEquationsToAdd = ReactionFormulas(contains(ReactionNames2,replace(rxnList(singleOrganismRxns{org}),'rn:','')));
        rxnNamesToAdd = singleStructs.(fieldNames{org}).rxns;
        for rxn=1:length(rxnEquationsToAdd)
            rxnEquationsToAdd(rxn) = replace(rxnEquationsToAdd(rxn),' +',['_',char(codeToAdd),' +']);
            rxnEquationsToAdd(rxn) = replace(rxnEquationsToAdd(rxn),' <',['_',char(codeToAdd),' <']);
            rxnEquationsToAdd(rxn) = cellstr([char(rxnEquationsToAdd(rxn)),'_',char(codeToAdd)]);
            rxnNamesToAdd(rxn) = cellstr([char(rxnNamesToAdd(rxn)),'_',char(codeToAdd)]);
        end
        communityRxnEquations(updateInterval) = rxnEquationsToAdd;
        communityRxnNames(updateInterval) = rxnNamesToAdd;
        counter = counter+length(updateInterval);
    end
    communityBounds = 20*ones(1,communityRxnNum);
    [throwAwayText,communityStruct] = evalc(char("createModel(communityRxnNames, communityRxnNames, communityRxnEquations', 'lowerBoundList', -communityBounds, 'upperBoundList', communityBounds);"));

    %Adding annotation fields to community reconstruction
    rxnCounter = 1;
    metCounter = 1;
    genCounter = 1;
    for org=1:nOrganisms
        %adding reaction-length annotations
        updateInterval = rxnCounter:(rxnCounter-1+length(singleStructs.(fieldNames{org}).rxns));
        communityStruct.grRules(updateInterval) = cellstr(singleStructs.(fieldNames{org}).grRules);
        communityStruct.rxnNames(updateInterval) = singleStructs.(fieldNames{org}).rxnNames;
        communityStruct.rxnECNumbers(updateInterval) = singleStructs.(fieldNames{org}).rxnECNumbers;
        communityStruct.rxnReferences(updateInterval) = singleStructs.(fieldNames{org}).rxnReferences;
        rxnCounter = rxnCounter+length(updateInterval);
        
        %adding metabolite-length annotations
        updateInterval = metCounter:(metCounter-1+length(singleStructs.(fieldNames{org}).mets));
        communityStruct.metNames(updateInterval) = singleStructs.(fieldNames{org}).metNames;
        metCounter = metCounter+length(updateInterval);
        
        %adding gene-length annotations
        updateInterval = genCounter:(genCounter-1+length(singleStructs.(fieldNames{org}).genes));
        communityStruct.genes(updateInterval) = cellstr(singleStructs.(fieldNames{org}).genes);
        communityStruct.geneNames(updateInterval) = singleStructs.(fieldNames{org}).geneNames;
        communityStruct.proteinNames(updateInterval) = singleStructs.(fieldNames{org}).proteinNames;
        communityStruct.proteins(updateInterval) = singleStructs.(fieldNames{org}).proteins;
        genCounter = genCounter+length(updateInterval);        
    end
    communityStruct.geneNames = communityStruct.geneNames';
    communityStruct.proteinNames = communityStruct.proteinNames';
    communityStruct.proteins = communityStruct.proteins';
    communityStruct.rxnECNumbers = communityStruct.rxnECNumbers';
    communityStruct.rxnReferences = communityStruct.rxnReferences';
    communityStruct = rmfield(communityStruct,'rxnGeneMat');
end

%% Identify components and give percentage of reactions in major component

fprintf("Identifying major components and finding percentage of reactions present in them.\n");

sMatNames = strings(0);
sMats = cell(0);
if ConsolidatedRecFlag
    sMatNames = [sMatNames "ConsolidatedRec"];
    sMats{length(sMatNames),1} = consolidatedStruct.S;
end
if SingleRecsFlag
    for org = 1:nOrganisms
        sMatNames = [sMatNames organismCodes(org)+" reconstruction"];
        sMats{length(sMatNames),1} = singleStructs.(fieldNames{org}).S;
    end
end
if CommunityRecFlag
    sMatNames = [sMatNames "CommunityRec"];
    sMats{length(sMatNames),1} = communityStruct.S;
end

for recIt=1:length(sMatNames)
    sMat = sMats{recIt,1};
    rxnsChecked = false(size(sMat,2),1);
    rxnsToCheckNext = false(size(sMat,2),1);
    rxnComponents = zeros(size(sMat,2),1);
    componentIndex = 0;
    for rxn = 1:size(sMat,2)
       nextRxn = rxn;
       if ~rxnsChecked(rxn)
           rxnsToCheckNext(rxn) = true;
           componentIndex = componentIndex+1;
           rxnComponents(rxn) = componentIndex;
           connectedCpds = find(sMat(:,nextRxn)~=0);
           connectedRxns = find(sum(abs(sMat(connectedCpds,:)))~=0);
           rxnsToCheckNext(connectedRxns) = true;
           rxnsToCheckNext(rxnsChecked) = false;
           while max(rxnsToCheckNext)==true
               rxnsToCheckNext(nextRxn) = false;
               rxnsChecked(nextRxn) = true;
               nextRxn = min(find(rxnsToCheckNext==1));
               connectedCpds = find(sMat(:,nextRxn)~=0);
               connectedRxns = find(sum(abs(sMat(connectedCpds,:)))~=0);
               rxnsToCheckNext(connectedRxns) = true;
               rxnsToCheckNext(rxnsChecked) = false;
               rxnComponents(connectedRxns) = componentIndex;
           end
       end
    end
    componentSizes = zeros(componentIndex,1);
    for component=1:componentIndex
       componentSizes(component) = length(find(rxnComponents==component));
    end
    networkData = struct();
    networkData.componentSizes = componentSizes;
    networkData.componentIndexes = rxnComponents;
    %printing:
    if sMatNames(recIt) == "CommunityRec"
        sortedComponentSizes = sort(componentSizes,'descend');
        fprintf("\tReactions in major components of the %-25s %5.2f%% \n", sMatNames(recIt)+":", 100*sum(sortedComponentSizes(1:nOrganisms))/sum(componentSizes));           
        disconnectedRxns = false(length(rxnComponents),1);
        disconnectedComponents = true(length(componentSizes),1);
        for i=1:length(componentSizes)
            if componentSizes(i)>=sortedComponentSizes(nOrganisms)
                disconnectedComponents(i) = false;
            end
        end
        for i=1:length(rxnComponents)
            if disconnectedComponents(rxnComponents(i))
                disconnectedRxns(i) = true;
            end
        end 
        networkData.disconnectedReactionIndexes = find(disconnectedRxns==true);
        networkData.percentageRxnsInMajorComponent = 100*sum(sortedComponentSizes(1:nOrganisms))/sum(componentSizes);
        networkData.numberRxnsInMajorComponent = sum(sortedComponentSizes(1:nOrganisms));
        communityStruct.networkData = networkData;
    else
        fprintf("\tReactions in major component of the %-26s %5.2f%% \n", sMatNames(recIt)+":", 100*max(componentSizes)/sum(componentSizes));
        disconnectedRxns = false(length(rxnComponents),1);
        disconnectedComponents = true(length(componentSizes),1);
        for i=1:length(componentSizes)
            if componentSizes(i)>=max(componentSizes)
                disconnectedComponents(i) = false;
            end
        end
        for i=1:length(rxnComponents)
            if disconnectedComponents(rxnComponents(i))
                disconnectedRxns(i) = true;
            end
        end
        networkData.disconnectedReactionIndexes = find(disconnectedRxns==true);
        networkData.percentageRxnsInMajorComponent = 100*max(componentSizes)/sum(componentSizes);
        networkData.numberRxnsInMajorComponent = max(componentSizes);
        if ConsolidatedRecFlag
            if recIt==1
                consolidatedStruct.networkData = networkData;
            else
                fieldNames = fieldnames(singleStructs);
                singleStructs.(fieldNames{recIt-1}).networkData = networkData;
            end
        else
            fieldNames = fieldnames(singleStructs);
            singleStructs.(fieldNames{recIt}).networkData = networkData;
        end
    end
end

%% Write out reconstruction(s) to file

if writeSBMLflag
    format = 'sbml';
    time = clock;
    timeString = strings(0);
    for t=1:length(time)-1
        if time(t)<10
            timeString = [timeString string(['0' char(string(time(t)))])];
        else
            timeString = [timeString string(char(string(time(t))))];
        end
    end
    timeString = [char(timeString(1)) '.' char(timeString(2)) '.' char(timeString(3)) '_' char(timeString(4)) '.' char(timeString(5))];
    printString = "Writing %s to file.\n";
    totalReactionsToPrint = 0;
    totalRecsToPrint = 0;
    if ConsolidatedRecFlag
        totalReactionsToPrint = totalReactionsToPrint+length(consolidatedStruct.rxns);
        totalRecsToPrint = totalRecsToPrint+1;
    end
    if SingleRecsFlag
        for org=1:nOrganisms
            totalReactionsToPrint = totalReactionsToPrint+length(singleStructs.(fieldNames{org}).rxns);
            totalRecsToPrint = totalRecsToPrint+1;
        end
    end
    if CommunityRecFlag
        totalReactionsToPrint = totalReactionsToPrint+length(communityStruct.rxns);
        totalRecsToPrint = totalRecsToPrint+1;
    end
    totalReactionsPrinted = 0;
    if totalRecsToPrint>1
        fprintf("Writing SBML files. This takes some time.\n");
        h = waitbar(0,'Printing first reconstruction...','Name','Writing SBML files');
    else
        fprintf("Writing SBML file. This takes some time.\n");
    end
    recsPrinted = 0;
	if ConsolidatedRecFlag
        recName = "ConsolidatedRec";
        if totalRecsToPrint>1
            waitbar(totalReactionsPrinted/totalReactionsToPrint,h,sprintf('Writing %s to file. %i of %i files written.',recName,recsPrinted,totalRecsToPrint));
        end
        fprintf(printString,recName);
        writeCbModel(consolidatedStruct,'format',format,'fileName',['ConsolidatedRec_' timeString]);
        totalReactionsPrinted = totalReactionsPrinted+length(consolidatedStruct.rxns);
        recsPrinted = recsPrinted+1;
    end
	if SingleRecsFlag
        for org=1:nOrganisms
            recName = string([char(organismCodes(org)), ' reconstruction']);
            if totalRecsToPrint>1
                waitbar(totalReactionsPrinted/totalReactionsToPrint,h,sprintf('Writing %s to file. %i of %i files written.',recName,recsPrinted,totalRecsToPrint));
            end
            fprintf(printString,recName);
            fieldNames = fieldnames(singleStructs);
            writeCbModel(singleStructs.(fieldNames{org}),'format',format,'fileName',[char(organismCodes(org)) 'Rec_' timeString]);
            totalReactionsPrinted = totalReactionsPrinted+length(singleStructs.(fieldNames{org}).rxns);
            recsPrinted = recsPrinted+1;
        end
    end
    if CommunityRecFlag
        recName = "CommunityRec";
        if totalRecsToPrint>1
            waitbar(totalReactionsPrinted/totalReactionsToPrint,h,sprintf('Writing %s to file. %i of %i files written.',recName,recsPrinted,totalRecsToPrint));
        end
        fprintf(printString,recName);
        writeCbModel(communityStruct,'format',format,'fileName',['CommunityRec_' timeString]);
        totalReactionsPrinted = totalReactionsPrinted+length(communityStruct.rxns);
        recsPrinted = recsPrinted+1;
    end
    if totalRecsToPrint>1
        close(h);
    end
end

outputStruct = struct();

if ConsolidatedRecFlag
    outputStruct.consolidatedStruct = consolidatedStruct;
end
if SingleRecsFlag
    outputStruct.singleStructs = singleStructs;
end
if CommunityRecFlag
    outputStruct.communityStruct = communityStruct;
end
if OrgRxnGenFlag
    for rxnnum=1:length(rxnOrganismGeneMatrix)
       reactionnumber = {};
       for orgnum=1:nOrganisms
           datafield = rxnOrganismGeneMatrix(rxnnum,orgnum);
           if char(datafield)>1
               if contains(datafield, ' | ')
                   datafield = strsplit(datafield);
                   numgenes=length(datafield(datafield~="|"));
                   reactionnumber{end+1} = numgenes;
               else
                   numgenes = length(datafield);
                   reactionnumber{end+1} = numgenes;
               end
           else
               numgenes = 0;
           end
       end
       rxnOrganismGeneMatrix(rxnnum,orgnum+1) = sum(~cellfun(@isempty,cellstr(rxnOrganismGeneMatrix(rxnnum,1:nOrganisms))));
       rxnOrganismGeneMatrix(rxnnum,orgnum+2) = sum(~cellfun(@isempty,cellstr(rxnOrganismGeneMatrix(rxnnum,1:nOrganisms))))/nOrganisms;
        if numgenes ~=0
           rxnOrganismGeneMatrix(rxnnum,orgnum+3) = strjoin(string(unique(cell2mat(reactionnumber))), ', ');
        else
           rxnOrganismGeneMatrix(rxnnum,orgnum+3) = 0;
        end
    end
    outputRxnOrgMat = cell(size(rxnOrganismGeneMatrix,1)+1,size(rxnOrganismGeneMatrix,2)+1);
    for i=1:nOrganisms
        outputRxnOrgMat{1,i+1} = cellstr(organismCodes(i));
    end
    outputRxnOrgMat{1,1} = cellstr("KEGG ID");
    outputRxnOrgMat{1,nOrganisms+2} = cellstr("Sum");
    outputRxnOrgMat{1,nOrganisms+3} = cellstr("Total");
    outputRxnOrgMat{1,nOrganisms+4} = cellstr("Genes");
    for i=1:length(rxnList)
        outputRxnOrgMat{i+1,1} = cellstr(replace(string(rxnList(i)),"rn:",""));
    end
    for i=1:size(rxnOrganismGeneMatrix,1)
       for j=1:size(rxnOrganismGeneMatrix,2)
           outputRxnOrgMat{i+1,j+1} = cellstr(rxnOrganismGeneMatrix(i,j));
       end
    end
    outputStruct.rxnOrganismGeneMatrix = outputRxnOrgMat;
end

if OmittedDataFlag
    omittedOutput = struct();
    omittedOutput.(char("OmittedReactions")) = omittedrxns;
    omittedOutput.(char("OmittedCompounds")) = omittedCPDs;
    outputStruct.omittedOutput = omittedOutput;
end

if DisconnectedReactionsFlag
    rows = 0;
    if ConsolidatedRecFlag
        rows = rows+1;
    end
    if SingleRecsFlag
        rows = rows+nOrganisms;
    end
    if CommunityRecFlag
        rows = rows+1;
    end
    disconnectedReactions = cell(rows,2);
    rowIt = 0;
    if ConsolidatedRecFlag
        rowIt = rowIt+1;
        disconnectedReactions(rowIt,1) = cellstr("ConsolidatedRec");
        disconnectedReactions(rowIt,2) = cellstr(strjoin(outputStruct.consolidatedStruct.rxns(outputStruct.consolidatedStruct.networkData.disconnectedReactionIndexes),"; "));
    end
    if SingleRecsFlag
        for i=1:nOrganisms
            rowIt = rowIt+1;
            disconnectedReactions(rowIt,1) = cellstr(organismCodes(i)+" reconstruction");
            fieldNames = fieldnames(singleStructs);
            disconnectedReactions(rowIt,2) = cellstr(strjoin(outputStruct.singleStructs.(fieldNames{i}).rxns(outputStruct.singleStructs.(fieldNames{i}).networkData.disconnectedReactionIndexes),"; "));
        end
    end
    if CommunityRecFlag
        rowIt = rowIt+1;
        disconnectedReactions(rowIt,1) = cellstr("CommunityRec");
        disconnectedReactions(rowIt,2) = cellstr(strjoin(outputStruct.communityStruct.rxns(outputStruct.communityStruct.networkData.disconnectedReactionIndexes),"; "));
    end
    outputStruct.disconnectedReactions = disconnectedReactions;
end

time = toc;

%% Final output notes
fprintf("\n\n\nAutoKEGGRec summary:\n\n")
fprintf("AutoKEGGRec completed its instructions in %.0f seconds. \n\n", time);
fprintf("%-55s %10s\n", "Organisms received as input:", string(nOrganisms));
if recBeingBuilt || OmittedDataFlag
    fprintf("%-55s %10d\n", "KEGG reactions assessed:", length(rxnsToBuildList));
    fprintf("\t%-51s %10d\n", "Reactions marked as redundant:", numGenReactions);
    fprintf("\t%-51s %10d\n", "Reactions marked as moved:", numMovedReactions);
    fprintf("\t%-51s %10d\n", "Reactions containing glycans:", sugarReactions  + nContainingReactions);
    fprintf("\t%-51s %10d\n", "Reactions involving polymerization:", polymerReactions);
    fprintf("\t%-51s %10d\n", "Reactions producing their own substrate:", reactionsProducingTheirOwnSubstrate);
    fprintf("%-55s %10d\n", "Unique compounds retrieved from KEGG:", length(cpdlist));
    fprintf("%-55s %10d\n", "Generic (massless) compounds identified:", length(genericCPDs));
    fprintf("%-55s %10d\n", "Reactions omitted for containing generic compounds:", length(reactionCheckList));
    fprintf("\n");
end
if recBeingBuilt
    if SingleRecsFlag
        fprintf("Single-organism reconstructions created for all the input organisms.\n");
    end
    if CommunityRecFlag
        fprintf("Community reconstruction created for all the input organisms.\n");
    end
    if ConsolidatedRecFlag
        fprintf("Consolidated/generic reconstruction created for all the input organisms.\n");
    end
    fprintf("\n");
end
if OrgRxnGenFlag
    fprintf("Organism-Reaction-Gene-Matrix added to the output structure.\n");
end
if OmittedDataFlag
    fprintf("Omitted reactions and compounds added to the output structure.\n");
end

fprintf("\nThanks for using AutoKEGGRec. Please cite our paper:\n E.Karlsen, C.Schulz and E.Almaas \n Automated generation of genome-scale metabolic draft reconstructions based on KEGG \n BMC Bioinformatics 2018 \n ");
end
