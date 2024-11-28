import fermi
import settings
import os
import pandas as pd


# test. change this to reading from file/table:


my_table = pd.read_csv(os.getenv("MASTERTABLE_PATH"), sep="\s+")


for source_name, ra, dec in zip(my_table["4FGLName"], my_table["RA"], my_table["DEC"]):
    source_name = source_name.replace("_", " ")

    my_source = fermi.Source(source_name, ra,  dec)

    my_source.get_fermi_lat_data(spacecraft=True)
    my_source.create_events_file()
    my_source.create_config_file()
