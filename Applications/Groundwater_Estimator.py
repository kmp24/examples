import os
import streamlit as st
from streamlit_folium import folium_static
import numpy as np
import datetime
from datetime import timedelta
import geopandas as gpd
import pandas as pd
from rasterstats import zonal_stats
import folium
from pydap.client import open_url
from pydap.cas.urs import setup_session
import xarray as xr

# Workaround for authentication issues
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# Function to handle PyDAP authentication
def openNetCDF(url, username, password):
    session = setup_session(username=username, password=password, check_url=url)
    store = xr.backends.PydapDataStore.open(url, session=session)
    ds = xr.open_dataset(store, decode_coords="all", mask_and_scale=False)
    return ds

# Function to get the file name based on the year, month, and day
def returnFileName(year, month, day):
    baseurl = 'https://hydro1.gesdisc.eosdis.nasa.gov/opendap/hyrax/GLDAS/GLDAS_CLSM025_DA1_D.2.2'
    filename = f'GLDAS_CLSM025_DA1_D.A{year}{month.zfill(2)}{day.zfill(2)}.022.nc4'
    urlbase = f"{baseurl}/{year}/{month.zfill(2)}/{filename}"
    return urlbase

# Function to calculate the difference between two dates
def calcDifference(earlierDate, laterDate):
    diff_calc = laterDate.GWS_tavg.data - earlierDate.GWS_tavg.data 
    diff_calc = diff_calc[0]  # Remove third dimension for visualization
    return diff_calc

# Function to get a folium basemap with ESRI tiles
def esriBasemap(miny, maxx):
    tile = "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}"
    attribu = 'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community'
    m = folium.Map([miny, maxx], zoom_start=2, tiles=tile, attr=attribu)
    return m

# Function to process uploaded file and save it temporarily
def save_uploaded_file(file_content, file_name):
    import tempfile
    import os
    import uuid

    _, file_extension = os.path.splitext(file_name)
    file_id = str(uuid.uuid4())
    file_path = os.path.join(tempfile.gettempdir(), "{}{}".format(file_id, file_extension))

    with open(file_path, "wb") as file:
        file.write(file_content.getbuffer())

    return file_path

# Function to convert date range to list of dates
def daterange(start, end):
    return [start + timedelta(n) for n in range(int((end - start).days))]

# Function to get bounding coordinates of a polygon
def boundsInPoly(df):
    minmax = ['minx', 'miny', 'maxx', 'maxy']
    for f in minmax:
        df[f] = df.bounds[f]
    return df[f]

def main():
    st.title('Groundwater Estimator Application')

    st.sidebar.markdown('''
    This application uses GLDAS Catchment Land Surface Model L4 daily 0.25 x 0.25 degree GRACE-DA1 V2.2 (GLDAS_CLSM025_DA1_D_2.2) as the
    basis for estimating change in groundwater between two dates.
    ''')

    option = st.sidebar.selectbox(
        'What area would you like to examine?',
        ('Global', 'Area of Interest')
    )

    st.write(option)

    if option == 'Area of Interest':
        uploaded_file = st.sidebar.file_uploader("Choose a file", type=["geojson", "gpkg"])

        file_path = save_uploaded_file(uploaded_file, uploaded_file.name)
        layer_name = os.path.splitext(uploaded_file.name)[0]

        infile = gpd.read_file(file_path)


    ### Sidebar for date selection
    st.sidebar.markdown('## Select Dates')
    st.sidebar.markdown('Select dates between Feb 1 2003 and 6 months prior to current date. The latency period for GLDAS/Grace can be up to 6 months.')
    # Supply the min and max dates available for GLDAS v2.2. 
    # There is a latency period of approximately 2-6 months, and data is updated mid-month.
    # (https://ldas.gsfc.nasa.gov/faq/gldas)
    # (https://hydro1.gesdisc.eosdis.nasa.gov/opendap/hyrax/GLDAS/GLDAS_CLSM025_DA1_D_EP.2.2/doc/README_GLDAS2.pdf)



    # min and max dates for each date input box
    today = datetime.date.today()
    min_date1, max_date1 = datetime.datetime(2004,2,2), datetime.datetime(today.year, today.month, today.day)
    min_date2, max_date2 = datetime.datetime(2004,2,2), datetime.datetime(today.year, today.month, today.day)

    prev_date = st.sidebar.date_input("Select prev date", min_value=min_date1, max_value=max_date1)
    later_date = st.sidebar.date_input("Select later date", min_value=min_date2, max_value=max_date2)

    st.sidebar.markdown('''The goal of the Global Land Data Assimilation System (GLDAS) is to ingest satellite- and ground-based observational data products, using advanced land surface modeling and data assimilation techniques, in order to generate optimal fields of land surface states and fluxes (Rodell et al., 2004a). https://ldas.gsfc.nasa.gov/gldas''')

    # Retrieve string year/month/date to select GLDAS files, pad single digit month and day with leading 0
    pyear, pmonth, pday = str(prev_date.year), str(prev_date.month).zfill(2), str(prev_date.day).zfill(2)
    lyear, lmonth, lday = str(later_date.year), str(later_date.month).zfill(2), str(later_date.day).zfill(2)

    st.echo()
    if st.button('Process Date'):
            with st.spinner("Data is processing"):
                # Calculate GWS change and add data to maps in the app
                if datetime.datetime(prev_date.year, prev_date.month, prev_date.day) != datetime.date.today:

                    if datetime.datetime(later_date.year, later_date.month, later_date.day) != datetime.date.today:

                        # Two dates for comparison
                        in_prev = returnFileName(pyear,pmonth,pday)
                        in_later = returnFileName(lyear, lmonth, lday)

                        # Open the netcdf files that correspond to the dates selected
                        ds1 = openNetCDF(in_prev, username, password)
                        ds2 = openNetCDF(in_later, username, password)

                        # get the lat/lon coordinate pairs in a format folium can use
                        lon, lat = np.meshgrid(ds1.lon.values.astype(np.float64), ds1.lat.values.astype(np.float64))

                        # Check the files exist
                        if xr.DataArray.notnull(ds1):

                            #st.write(in_prev)

                            if xr.DataArray.notnull(ds2):
                                
                                # Calculate the difference between the two files
                                changeGWS = calcDifference(ds1, ds2)

                                # Visualize the change using folium
                                img = np.flipud(changeGWS)
                                img = (img/10).astype(int) # convert mm to cm and round for visualization
                                img = np.ma.masked_where(img==0, img)

                                import matplotlib
                                cmap = matplotlib.cm.seismic
                                cmap.set_bad('grey', alpha=0)

                                if option == 'Global':
                                    ## Folium map
                                    # Initialize a basemap for viewing the comparison
                                    esriBasemap(15.0, 0.0)
                                    folium.raster_layers.ImageOverlay(img, mercator_project=True, colormap=cmap, opacity=.5,
                                                                bounds =[[lat.min(), lon.min()], [lat.max(), lon.max()]]).add_to(m)

                                    folium_static(m)


                                elif option == 'Area of Interest':
                                    
                                    data = pd.DataFrame()
                                    array = img
                                    affine = rasterio.transform.from_origin(lon.min(), lat.max(), 0.25, 0.25)

                                    # Look at nodata/0

                                    # Zonal stats method -  consider all_touched=True/False
                                    GWS_stats = zonal_stats(infile, array, affine=affine, stats=['count','sum', 'mean'])
                                    gwsstats = pd.DataFrame(GWS_stats, columns={'count','sum','mean'})

                                    # add to dataframe
                                    data = data.append(gwsstats) # do something with these

                                    # Initialize a basemap for viewing the comparison
                                    esriBasemap(infile.geometry.centroid.y.iloc[0],infile.geometry.centroid.x.iloc[0])
                                    
                                    ## Folium map
                                    folium.raster_layers.ImageOverlay(img, mercator_project=True, colormap=cmap, opacity=.5,
                                                                bounds=[[lat.min(), lon.min()], [lat.max(), lon.max()]]).add_to(m)

                                    # Add AOI polygon
                                    for _, r in infile.iterrows():
                                        # From Geopandas documentation
                                        sim_geo = gpd.GeoSeries(r['geometry']).simplify(tolerance=0.001)
                                        geo_j = sim_geo.to_json()
                                        geo_j = folium.GeoJson(data=geo_j)
                                        geo_j.add_to(m)

                                    folium_static(m)
                            else:
                                st.write('Data is empty, check dates')    

                        else:
                            st.write('Data is empty, check dates')           

                    else:
                        st.write('Dates selected are not available due to latency of the GLDAS product')

                else:
                        st.write('Dates selected are not available due to latency of the GLDAS product')


if __name__ == "__main__":
    main()