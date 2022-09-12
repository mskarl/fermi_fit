import fermi

# test. change this to reading from file/table:

my_source = fermi.Source('4FGL J0103.5+1526', 15.8584, 15.4402)
# my_source = fermi.Source('4FGLJ0103.5+1526', 15.6, 15.6)

# # get data. change this to a parsing option maybe...

my_source.get_fermi_lat_data(spacecraft=True)
# #my_source.get_spacecraft_data()
my_source.create_events_file()
my_source.create_config_file()
my_source.setup_analysis()
my_source.free_parameter()
my_source.fit_llh()
my_source.calculate_bowtie()
my_source.calculate_sed()
