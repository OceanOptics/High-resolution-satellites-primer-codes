function data = importfile_SPMfromRWS(filename, dataLines)
%IMPORTFILE Import data from a text file
%  UNTITLED = IMPORTFILE(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  UNTITLED = IMPORTFILE(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  Untitled = importfile("/path/to/data/SPM_fromRWS/data.txt", [2, Inf]); %data.txt is data downloaded directly from Rijkwaterstaad (waterinfo) database
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 14-Feb-2022 13:18:08

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 50);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["MONSTER_IDENTIFICATIE", "MEETPUNT_IDENTIFICATIE", "LOCATIE_CODE", "TYPERING_OMSCHRIJVING", "TYPERING_CODE", "GROOTHEID_OMSCHRIJVING", "GROOTHEID_CODE", "PARAMETER_OMSCHRIJVING", "PARAMETER_CODE", "CAS_NR", "EENHEID_CODE", "HOEDANIGHEID_OMSCHRIJVING", "HOEDANIGHEID_CODE", "COMPARTIMENT_OMSCHRIJVING", "COMPARTIMENT_CODE", "WAARDEBEWERKINGSMETHODE_OMSCHRIJVING", "WAARDEBEWERKINGSMETHODE_CODE", "WAARDEBEPALINGSMETHODE_OMSCHRIJVING", "WAARDEBEPALINGSMETHODE_CODE", "BEMONSTERINGSSOORT_OMSCHRIJVING", "BEMONSTERINGSSOORT_CODE", "WAARNEMINGDATUM", "WAARNEMINGTIJDMETCET", "LIMIETSYMBOOL", "NUMERIEKEWAARDE", "ALFANUMERIEKEWAARDE", "KWALITEITSOORDEEL_CODE", "REFERENTIE", "NOTITIE_CODE", "NOTITIE_OMSCHRIJVING", "STATUSWAARDE", "OPDRACHTGEVENDE_INSTANTIE", "MEETAPPARAAT_OMSCHRIJVING", "MEETAPPARAAT_CODE", "BEMONSTERINGSAPPARAAT_OMSCHRIJVING", "BEMONSTERINGSAPPARAAT_CODE", "PLAATSBEPALINGSAPPARAAT_OMSCHRIJVING", "PLAATSBEPALINGSAPPARAAT_CODE", "BEMONSTERINGSHOOGTE", "REFERENTIEVLAK", "EPSG", "X", "Y", "ORGAAN_OMSCHRIJVING", "ORGAAN_CODE", "TAXON_NAME", "GROEPERING_OMSCHRIJVING", "GROEPERING_CODE", "GROEPERING_KANAAL", "GROEPERING_TYPE"];
opts.VariableTypes = ["string", "categorical", "categorical", "string", "string", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "string", "string", "categorical", "categorical", "string", "string", "categorical", "double", "categorical", "categorical", "datetime", "datetime", "string", "double", "string", "categorical", "string", "string", "string", "categorical", "categorical", "string", "string", "categorical", "double", "string", "string", "double", "categorical", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["MONSTER_IDENTIFICATIE", "TYPERING_OMSCHRIJVING", "TYPERING_CODE", "HOEDANIGHEID_OMSCHRIJVING", "HOEDANIGHEID_CODE", "WAARDEBEWERKINGSMETHODE_OMSCHRIJVING", "WAARDEBEWERKINGSMETHODE_CODE", "LIMIETSYMBOOL", "ALFANUMERIEKEWAARDE", "REFERENTIE", "NOTITIE_CODE", "NOTITIE_OMSCHRIJVING", "MEETAPPARAAT_OMSCHRIJVING", "MEETAPPARAAT_CODE", "PLAATSBEPALINGSAPPARAAT_OMSCHRIJVING", "PLAATSBEPALINGSAPPARAAT_CODE", "ORGAAN_OMSCHRIJVING", "ORGAAN_CODE", "TAXON_NAME", "GROEPERING_OMSCHRIJVING", "GROEPERING_CODE", "GROEPERING_KANAAL", "GROEPERING_TYPE"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["MONSTER_IDENTIFICATIE", "MEETPUNT_IDENTIFICATIE", "LOCATIE_CODE", "TYPERING_OMSCHRIJVING", "TYPERING_CODE", "GROOTHEID_OMSCHRIJVING", "GROOTHEID_CODE", "PARAMETER_OMSCHRIJVING", "PARAMETER_CODE", "CAS_NR", "EENHEID_CODE", "HOEDANIGHEID_OMSCHRIJVING", "HOEDANIGHEID_CODE", "COMPARTIMENT_OMSCHRIJVING", "COMPARTIMENT_CODE", "WAARDEBEWERKINGSMETHODE_OMSCHRIJVING", "WAARDEBEWERKINGSMETHODE_CODE", "WAARDEBEPALINGSMETHODE_OMSCHRIJVING", "BEMONSTERINGSSOORT_OMSCHRIJVING", "BEMONSTERINGSSOORT_CODE", "LIMIETSYMBOOL", "ALFANUMERIEKEWAARDE", "KWALITEITSOORDEEL_CODE", "REFERENTIE", "NOTITIE_CODE", "NOTITIE_OMSCHRIJVING", "STATUSWAARDE", "OPDRACHTGEVENDE_INSTANTIE", "MEETAPPARAAT_OMSCHRIJVING", "MEETAPPARAAT_CODE", "BEMONSTERINGSAPPARAAT_OMSCHRIJVING", "PLAATSBEPALINGSAPPARAAT_OMSCHRIJVING", "PLAATSBEPALINGSAPPARAAT_CODE", "REFERENTIEVLAK", "ORGAAN_OMSCHRIJVING", "ORGAAN_CODE", "TAXON_NAME", "GROEPERING_OMSCHRIJVING", "GROEPERING_CODE", "GROEPERING_KANAAL", "GROEPERING_TYPE"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "WAARNEMINGDATUM", "InputFormat", "dd/MM/yyyy");
opts = setvaropts(opts, "WAARNEMINGTIJDMETCET", "InputFormat", "HH:mm:ss");
opts = setvaropts(opts, "WAARDEBEPALINGSMETHODE_CODE", "TrimNonNumeric", true);
opts = setvaropts(opts, "WAARDEBEPALINGSMETHODE_CODE", "ThousandsSeparator", ",");

% Import the data
data = readtable(filename, opts);

end
