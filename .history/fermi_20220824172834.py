import sys
import os
import settings
import time
import yaml
import numpy as np
import glob

import requests
from mechanize import Browser
from astropy.io import fits
from scipy.interpolate import interp1d

from fermipy.gtanalysis import GTAnalysis


class Source:

    def __init__(self, source_name, ra, dec, alert_time=None):
        print(source_name, ra, dec, alert_time, os.getenv("OUTPUT_DIR"))
        self.source_name = source_name
        self.ra = ra
        self.dec = dec
        self.alert_time = alert_time
        self.working_dir = os.path.join(os.getenv("OUTPUT_DIR"), self.source_name.replace(" ", ""))
        self.data_dir = os.path.join(self.working_dir, "fermi_data")
        self.spacecraft_dir = os.path.join(os.getenv("OUTPUT_DIR"), "spacecraft")
        print("output will be saved to ", self.working_dir)

        self.gta = None
        self.free_params = False
        return 
        

    # TODO implement time stuff.. or maybe apply cut later?
    def get_fermi_lat_data(self, spacecraft=True): # from Theo's skript



        def get_download_links(html):
            split = html.decode().split('wget')
            status = int(html.decode().split('he state of your query is ')[1][:1])
            if status == 2:
                return [(i.split('</pre>'))[0].strip().replace('\n', '')
                        for i in split[2:]]
            else:
                return []

        def download_file(url, outfolder):
            local_filename = os.path.join(outfolder, url.split('/')[-1])
            # NOTE the stream=True parameter below
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(local_filename, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192): 
                        # If you have chunk encoded response uncomment if
                        # and set chunk_size parameter to None.
                        #if chunk: 
                        f.write(chunk)
            return

        url = "https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi"

        br = Browser()
        br.set_handle_robots(False)
        br.open(url)
        br.select_form(nr=1)
        br["coordfield"] = str(self.ra) + ", " + str(self.dec)
        br["coordsystem"] = [u'J2000']
        br["timefield"] = "START, END"
        br["shapefield"] = "15"
        br["energyfield"] = "100, 1000000"
        # radius br["shapefield"]
        # time br["timefield"], br["timetype"]
        # we load the spacecraft file separately
        # br.form.find_control('spacecraft').items[0].selected = False

        response = br.submit()
        r_text = response.get_data()
        query_url = r_text.decode().split('at <a href="')[1].split('"')[0]

        print('Query URL {}'.format(query_url))
        seconds = float(r_text.decode().split('complete is ')[1].split(' ')[0])
        wait = 0.75 * seconds
        print('Wait at least {} seconds for the query to complete'.format(wait))
        time.sleep(wait)

        html = requests.get(query_url).text.encode('utf8')
        download_urls = get_download_links(html)
        while len(download_urls) == 0:
            print('Query not yet finished...Wait 10 more seconds')
            time.sleep(10)
            html = requests.get(query_url).text.encode('utf8')
            download_urls = get_download_links(html)

        if not os.path.exists(self.data_dir):
            os.makedirs(self.data_dir)
        
        # remove old data

        old_files = glob.glob(os.path.join(self.data_dir, '*.fits'))
        for f in old_files:
            print("delete ", f)
            os.remove(f)

        for tmp_url in download_urls:
            download_file(tmp_url, self.data_dir)

        return


    def get_spacecraft_data(self):
        os.system("wget -N -P " + self.spacecraft_dir + \
             " https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/mission/spacecraft/lat_spacecraft_merged.fits")


    def create_events_file(self):
        
        files = os.listdir(self.data_dir)
        ph_files = [os.path.join(self.data_dir, f) for f in files if "PH" in f]
        with open(os.path.join(self.data_dir, "events.txt"), 'w+') as event_file:
            event_file.write("\n".join(ph_files))
        return


    # def get_recent_catalog(self):
    #     os.system("wget -N -P " + os.path.join(os.getenv("OUTPUT_DIR"), "catalogs") + " )
    #     return


    # def get_latest_diffuse_models(self):  diffuse models are included in the fermi package stuff. this 
    # might need to be updated. 
    #     return 


#make_cuts()


#create_xml_file

    def get_time_window(self):
        files = os.listdir(self.data_dir)
        files = [f for f in files if 'PH' in f]
        tmin = []
        tmax = []
        for f in files:
            x = fits.open(os.path.join(self.data_dir, f))
            tmin.append(x[1].header['TSTART'])
            tmax.append(x[1].header['TSTOP'])
        return float(np.min(tmin)), float(np.max(tmax))


    def create_config_file(self):
        config_path = os.path.join(self.data_dir, "config.yaml")

        with open(os.path.join(".", "default.yaml"), "r") as stream:
            config = yaml.safe_load(stream)

        config['selection']['target'] = self.source_name

        config['selection']['ra'] = self.ra
        config['selection']['dec'] = self.dec
        config['data']['evfile'] = os.path.join(self.data_dir, "events.txt")
        spacecraft_dir = glob.glob(os.path.join(self.data_dir, "*SC*.fits"))
        config['data']['scfile'] = spacecraft_dir[0]

        # replace latest catalog here in case
        config['model']['catalogs'] = [os.path.join(os.getenv("OUTPUT_DIR"), "catalogs/gll_psc_v30.fit")]

        tmin, tmax = self.get_time_window()

        config['selection']['tmin'] = tmin
        config['selection']['tmax'] = tmax

        with open(config_path, 'w+') as stream:
            config = yaml.dump(config, stream, default_flow_style=False)



    def setup_analysis(self, config_file="config.yaml"):

        # setup the analysis
        config_path = os.path.join(self.data_dir, config_file)
        print("initialize analysis") # TODO change prints to logger?
        self.gta = GTAnalysis(config_path, logging={'verbosity': 3})
        print("setup analysis")
        self.gta.setup()
        print("optimize")
        self.gta.optimize()

        self.free_params = False
        return


    def free_parameter(self):
        # free parameter 
        if not self.free_params:
            print("free parameter")
            # TODO are these the correct steps? ask paolo
            self.gta.free_sources(distance=3.0, minmax_ts=[2, None], exclude=['isodiff', 'galdiff'], 
            pars="norm")
            self.gta.free_source("galdiff")
            self.gta.free_source("isodiff")
            if self.source_name:
                self.gta.free_source(self.source_name)

        self.free_params == True
        return 


    def fit_llh(self, filename="llh.npy"):
        self.free_parameter()
        print("fit likelihood ", self.working_dir + filename)
        self.gta.fit(retries=50) # TODO theo has a retries argument here
        self.gta.write_roi(os.path.join(self.working_dir, filename))
        return

    
    def calculate_bowtie(self, filename="bowtie.npy"):
        print("calculate bowtie for " + self.source_name, os.path.join(self.working_dir, filename))
        bowtie = self.gta.bowtie(self.source_name)
        np.save(os.path.join(self.working_dir, filename), bowtie)
        return


    def calculate_sed(self, filename="sed.fits"):
        print("calculate the sed ", os.path.join(self.working_dir, filename))
        self.gta.sed(self.source_name, outfile=os.path.join(self.working_dir, filename), 
                    free_radius=self.__get_95_psf(100), free_background=True, 
                    loge_bins=list(np.arange(np.log10(100),
                                    np.log10(1000000), 0.5)))
        # theo has some more parameter: 
        # loge_bins=list(np.arange(np.log10(args['emin']),
        #                            np.log10(args['emax']), 0.5)),
        #             free_pars=free_pars,
        #             free_radius=args['free_radius'],
        #             free_background=True
        return


    def __get_68_psf(self, E):
        x = np.genfromtxt(os.path.join(os.getenv("PSF_DIR"), 'lat_68_psf.txt'), delimiter = ',')
        return interp1d(x[:,0], x[:,1])(E)    
    
    def __get_95_psf(self, E):
        x = np.genfromtxt(os.path.join(os.getenv("PSF_DIR"), 'lat_95_psf.txt'), delimiter = ',')
        return interp1d(x[:,0], x[:,1])(E)