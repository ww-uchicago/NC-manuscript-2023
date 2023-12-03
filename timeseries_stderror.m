function sem = timeseries_stderror(data)

sem = std(data,0,2)/sqrt(size(data,2));
end