def plot_max_point_path(stm, etm, read_si=True, read_wic=False,
                        mlt_range=[8, 16], mlat_range=[70, 85]):

    from track_max_point import find_filenames, read_wic_si_data
    from track_max_point import find_max_intensity_point
    from matplotlib.dates import DateFormatter

    dtms, fnames = find_filenames(stm, etm, read_si=read_si, read_wic=read_wic)
    mlats = []
    mlts = []
    vals = [] 
    for dtm in dtms:            
        fname, data = read_wic_si_data(dtm, read_si=read_si, read_wic=read_wic)
        (maxval, maxloc) = find_max_intensity_point(data, mlt_range=mlt_range,
                                                    mlat_range=mlat_range)
        mlats.append(data['mlat'][maxloc]) 
        mlts.append(data['mlt'][maxloc])
        vals.append(maxval)
            
    fig, ax = plt.subplots()
    ax.plot_date(dtms, mlts)

    return fig

# run the code
import datetime as dt
import matplotlib.pyplot as plt

#read_si=True; read_wic=False
read_si=False; read_wic=True
mlt_range = [9, 15]
mlat_range = [70, 85]
dtm = dt.datetime(2002, 3, 18, 16, 02)
stm = dt.datetime(2002, 3, 18, 15, 00)
etm = dt.datetime(2002, 3, 18, 17, 22)
fig = plot_max_point_path(stm, etm, read_si=True, read_wic=False,
                          mlt_range=[8, 16], mlat_range=[70, 85])
plt.show()

