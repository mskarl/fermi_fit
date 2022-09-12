import fermi
import settings
import os
import pandas as pd


# test. change this to reading from file/table:


my_table = pd.read_csv(os.getenv("MASTERTABLE_PATH"), sep="\s+")

for item in zip(my_table["4FGLName"].values, my_table["RA"].values, my_table["DEC"].values):
    print(item)


for source_name, ra, dec in zip(my_table["4FGLName"], my_table["RA"], my_table["DEC"]):
    source_name = source_name.replace("_", " ")
    print(source_name, ra, dec)

    my_source = fermi.Source(source_name, ra,  dec)

# get data. change this to a parsing option maybe...
    if False:
        my_source.get_fermi_lat_data(spacecraft=True)
        my_source.create_events_file()
        my_source.create_config_file()
        # my_source.setup_analysis()
        # my_source.free_parameter()
        # my_source.fit_llh()
        # my_source.calculate_bowtie()
        # my_source.calculate_sed()
