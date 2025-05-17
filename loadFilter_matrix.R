library(data.table)

dt = fread('Resources/globalterrorismdb.csv')
names(dt)
dt[1:5, 'country']
dt = dt[country == 217]
fwrite(dt, 'Resources/countryFil_globalterrorismdb.csv')
