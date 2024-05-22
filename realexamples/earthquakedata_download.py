import pandas as pd
import numpy as np
import urllib.request
import shutil
import os

def download_data(year_start, year_end):
    url = ("https://earthquake.usgs.gov/fdsnws/event/1/query.csv?"
           f"starttime={year_start}-01-01%2000:00:00&endtime={year_end}"
           "-12-31%2023:59:59&maxlatitude=46&minlatitude=22&maxlongitude=150"
           "&minlongitude=122&minmagnitude=2.5&eventtype=earthquake&orderby=time-asc")

    filename = f"earthquakes/{year_start}_{year_end}.csv"
    os.makedirs("earthquakes", exist_ok=True)
    urllib.request.urlretrieve(url, filename)
    return filename

if __name__ == "__main__":
    filenames = []
    filenames.append(download_data(1990, 1999))
    filenames.append(download_data(2000, 2009))
    filenames.append(download_data(2010, 2019))

    # combine csv files
    with open('earthquakes/1990_2019.csv', 'wb') as wfd:
        for i, f in enumerate(filenames):
            with open(f, 'rb') as fd:
                if i != 0:
                    fd.readline()  # Throw away header on all but first file
                shutil.copyfileobj(fd, wfd)

df = pd.read_csv("earthquakes/1990_2019.csv")
# filter earthquakes in Japan
df = df[df["place"].str.contains("Japan")]
df = df[["time", "longitude", "latitude"]]
df2=df[ df['time']>'2019-01-01 00:00:00']
# write time as a days since 2019-01-01
basedate = pd.Timestamp(("2019-01-01T00:00:00Z"))
df2["time"] = pd.to_datetime(df2["time"])
df2["time"] = df2["time"].apply(lambda x: (x - basedate).total_seconds() / 60 / 60 / 24)  # days
t=np.array(df2['time'])
np.savez('earthquakes/earthquakes.npz', a=t)
