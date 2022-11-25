function [SPM_table,stations,data, Lat, Lon] = process_fieldSPM(filename)

data = importfile_SPMfromrws(filename);

SPM = data.NUMERIEKEWAARDE(2:end);
SPM(SPM >= 999) = NaN;

data.WAARNEMINGTIJDMETCET.Format = 'dd.MM.uuuu HH:mm';
data.WAARNEMINGDATUM.Format = 'dd.MM.uuuu';
date = data.WAARNEMINGDATUM; date= date(2:end);
time = data.WAARNEMINGTIJDMETCET; time = time(2:end);
time.Format = 'HH:mm:ss'; %time = char(time);
date_table = table(date);

%%coordinates
[Lat, Lon] = utm2deg(data.X(2:end), data.Y(2:end), repmat('31 N',length(data.X(2:end)),1));

SPMtable = date + timeofday(data.WAARNEMINGTIJDMETCET(2:end));
SPMtable.Format = 'dd.MM.uuuu';

myData = array2table(SPM);
myData{:,2} = SPMtable(:);
myData{:,3} = Lat(:);
myData{:,4} = Lon(:);
myData{:,5} = time(:);
myData{:,6} = data.MEETPUNT_IDENTIFICATIE(2:end);

myData.Properties.VariableNames = {'SPM','datetime','Lat','Lon','time','Stations'};
SPM_table = table2timetable(myData);

SPM_table.datetime.Format = 'dd.MM.uuuu';

datastatus = data(2:end,5);
stations = categories(data.MEETPUNT_IDENTIFICATIE);

end

