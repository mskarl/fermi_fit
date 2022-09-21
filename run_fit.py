import fermi
import settings
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--source", default=0, type=int, help="index in master table")
parser.add_argument("--data_dir", default=None, type=str)
parser.add_argument("--output_dir", default=None, type=str)
parser.add_argument("--force_power_law", default=False, type=bool, action="store_true", help="Sometimes the spectral " + \
                    "type in the catalog is a log parabola. This forces the model to use a power law.")

args = parser.parse_args()
print("my arguments")
print(args)

my_table = pd.read_csv(os.getenv("MASTERTABLE_PATH"), sep="\s+")

source_name, ra, dec = my_table["4FGLName"][args.source], my_table["RA"][args.source], my_table["DEC"][args.source]
source_name = source_name.replace("_", " ")

my_source = fermi.Source(source_name, ra,  dec, data_dir=args.data_dir, output_dir=args.output_dir)

# get data. change this to a parsing option maybe...
# my_source.get_fermi_lat_data(spacecraft=True)
# my_source.create_events_file()
# my_source.create_config_file()
if args.force_power_law:
    my_source.add_power_law_source()
    
my_source.setup_analysis()
my_source.free_parameter()
my_source.fit_llh()
my_source.calculate_bowtie()
my_source.calculate_sed()