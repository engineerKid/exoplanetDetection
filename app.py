import os
import logging
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive
from astroquery.mast import Catalogs
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np                              # used for numerical operations.
from datetime import datetime
import lightkurve as lk    # analyzing light curve data from missions like TESS and Kepler.

# packages for frontend
from flask import Flask, jsonify, render_template, request


#######################################
#  Constants used in the applicaiton
########################################

# Base directory for data storage
BASE_DATA_DIR = os.path.join(os.getcwd(), 'data')
os.makedirs(BASE_DATA_DIR, exist_ok=True)

# Define sky regions
RA_CENTER = 290  # in degrees
DEC_CENTER = 44  # in degrees
SEARCH_RADIUS = 0.2 # in degrees 

# Program assumes that stars/exoplanets are related if within this distance
SEPARATION_ARCSEC = 1 

# Define stellar filter ranges as constants
TMAG_MIN = 4
TMAG_MAX = 16
RAD_MIN = 0.5  # Minimum stellar radius
RAD_MAX = 10  # Maximum stellar radius
TEFF_MIN = 3000  # Minimum effective temperature
TEFF_MAX = 8000  # Maximum effective temperature, 8k good limit


##########################################
# START: SET API TO SERVE FRONTEND CALLS
##########################################

app = Flask(__name__)

# Default API endpoint
@app.route('/')
def index():
    # return render_template('index.html')

    # passing coordinates to frontend for Plotly.js to lock panning/zooming beyound these boundaries 
    # as there is no data beyond that
    return render_template('index.html', ra_center=RA_CENTER, dec_center=DEC_CENTER, search_radius=SEARCH_RADIUS)


# API endpoint for getting Star/Exoplanet data
@app.route('/get_star_data')
def get_star_data():

    # Call prepare_dataset directly, which ensures the file is created if it doesn't exist
    star_df = prepare_dataset(RA_CENTER, DEC_CENTER, SEARCH_RADIUS)

    # Check if star_df is None, which indicates that prepare_dataset failed
    if star_df is None:
        return jsonify({"error": "Failed to prepare dataset"}), 500

   # Replace NaN values with None so they can be properly serialized into JSON
    star_df = star_df.replace({np.nan: None})
   
    # We should update the return for a selection of specific columns instead of entire dataset, e.g. below
    # return jsonify(star_df[['ra', 'dec', 'Tmag', 'c_has_exoplanet', 'Teff', 'lum']].to_dict(orient='records'))  

    # Converts the entire DataFrame (star_df) to a dictionary and returns all the columns as JSON.
    return jsonify(star_df.to_dict(orient='records'))   
    

# API endpoint for folded light curve and time-series light curve data
@app.route('/analyze_star', methods=['POST'])
def analyze_star():
    tic_id = request.json.get('tic_id')
    if not tic_id:
        return jsonify({'error': 'No TIC ID provided.'}), 400

    # Analyze the light curve
    result = analyze_light_curve(tic_id)
    if result is None:
        return jsonify({'error': f'Analysis failed for TIC ID {tic_id}.'}), 500
    else:
        return result


##########################################
# END: SET API TO SERVE FRONTEND CALLS
##########################################


# Configure logging
logging.basicConfig(level=logging.INFO)

####################################################################################
# classify_star
# Apply temperature ranges to assign spectral types (O, B, A, F, G, K, M).
####################################################################################
 
def classify_star(teff):
    # Pandas library function: 
    #   returns True if teff is "not available" (NA) or NaN (Not a Number)
    if pd.isna(teff):
        return 'Unknown'

    if teff > 30000:
        return 'O-type'
    elif teff > 10000:
        return 'B-type'
    elif teff > 7500:
        return 'A-type'
    elif teff > 6000:
        return 'F-type'
    elif teff > 5200:
        return 'G-type'
    elif teff > 3700:
        return 'K-type'
    else:
        return 'M-type'
        

####################################################################################
# cross_match_coordinates
# Cross-match exoplanet and star coordinates within the given separation
####################################################################################

def cross_match_coordinates(exoplanet_df, star_df, separation_arcsec=SEPARATION_ARCSEC):
    # Ensure valid RA and DEC for SkyCoord
    exoplanet_df.dropna(subset=['ra', 'dec'], inplace=True)
    star_df.dropna(subset=['ra', 'dec'], inplace=True)

    # Convert to SkyCoord objects with proper units
    exoplanet_coords = SkyCoord(ra=exoplanet_df['ra'].values * u.deg, dec=exoplanet_df['dec'].values * u.deg)
    star_coords = SkyCoord(ra=star_df['ra'].values * u.deg, dec=star_df['dec'].values * u.deg)

    # Cross-match within a given separation (default 1 arcsecond)
    idx, d2d, _ = star_coords.match_to_catalog_sky(exoplanet_coords)
    max_sep = separation_arcsec * u.arcsec
    sep_constraint = d2d < max_sep

    return idx, sep_constraint, exoplanet_df


#########################################################
# # Observation Details - Information for UI info box 
########################################################
def meta_observation_details(lc):
    # Extract info from lc - striclty for showing data in the UI info box - put in JSON format
    start_bjd = round(lc.time.bkjd[0])
    end_bjd = round(lc.time.bkjd[-1])
    start_date = datetime.fromisoformat(lc.time[0].isot).strftime("%b %d, %Y")  # ISO format - Result: "Nov 15, 2018"
    end_date = datetime.fromisoformat(lc.time[-1].isot).strftime("%b %d, %Y")  # Result: "Nov 15, 2018"
    observation_duration = int(lc.time[-1].value - lc.time[0].value)  # In days, just get the integer value
    cadence_interval = np.median(np.diff(lc.time.value))  # In days
    cadence_interval_minutes = round(cadence_interval * 24 * 60) # in minutues

    # JSON - Observation Details (Data file info)
    meta_observation_details = {
        'source': "TESS / SPOC",
        'observation_start_bjd': start_bjd, # (start date in days since TESS mission start)
        'observation_end_bjd': end_bjd, # (end date in days since TESS mission start)""
        'observation_start_date': start_date,
        'observation_end_date': end_date,
        'observation_duration': observation_duration,
        'cadence_interval_minutes': cadence_interval_minutes,
    }
    # logging.info("Data File Information: %s", meta_observation_details)
    return meta_observation_details


#########################################################
# # Stellar Information: Information for UI info box 
#########################################################
def meta_stellar_details(tic_id):
    csv_file = f"data/star_tic_{tic_id}.csv"  # Adjust the path and filename format as needed

    # Load the CSV file into a DataFrame
    star_df = pd.read_csv(csv_file)
    logging.info(f"Listing all the columns for meta_stellar_details: {(star_df.columns)}")
    logging.info(f"Listing row count for meta_stellar_details: {len(star_df)}")

    # logging.info("All columns in the dataset:")
    # for column in star_df.columns:
    #     logging.info(column)

    # Define the key columns you want to extract
    key_columns = ['ra', 'dec', 'Tmag', 'Teff', 'rad', 'mass', 'lum', 'logg', 'd', 'pmRA', 'pmDEC']

    # Filter to only include columns that are present in the DataFrame
    available_columns = [col for col in key_columns if col in star_df.columns]
    
    # Drop rows with any missing values in the specified columns and select the first valid row
    filtered_df = star_df.dropna(subset=available_columns)
    logging.info(f"Listing row count for filtered_df: {len(filtered_df)}")

    # Return the first row as a dictionary if available, or an error message if not
    if not filtered_df.empty:
        return filtered_df[available_columns].iloc[0].to_dict()

    logging.warning(f"No valid data available for TIC ID {tic_id}")
    return {"error": "No valid data available"}    


####################################################################################
# download_exoplanet_data
# Download exoplanet data from NASA Exoplanet Archive
####################################################################################

def download_exoplanet_data(ra_center, dec_center, search_radius):
    """Download exoplanet data for the specified region and save locally."""
    exoplanet_filename = f"exoplanet_ra{ra_center}_dec{dec_center}_radius{search_radius}.csv"
    exoplanet_file = os.path.join(BASE_DATA_DIR, exoplanet_filename)

    if os.path.exists(exoplanet_file):
        logging.info(f"Exoplanet data for the region already exists locally.")
        exoplanet_df = pd.read_csv(exoplanet_file)
    else:
        logging.info(f"Downloading exoplanet data for RA={ra_center}, Dec={dec_center}, radius={search_radius} deg from NASA Exoplanet Archive.")
        coord = SkyCoord(ra=ra_center * u.deg, dec=dec_center * u.deg)
        radius = search_radius * u.deg

        try:
            exoplanet_table = NasaExoplanetArchive.query_region(
                coordinates=coord,
                radius=radius,
                table="pscomppars"
            )
        except Exception as e:
            logging.error(f"An error occurred during NASA Exoplanet Archive Catalog query: {e}")
            return None

        exoplanet_df = exoplanet_table.to_pandas()
        logging.info(f"Number of exoplanets found: {len(exoplanet_df)}")

        # Check if 'sky_coord.ra' and 'sky_coord.dec' exist and have valid values
        if 'sky_coord.ra' in exoplanet_df.columns and 'sky_coord.dec' in exoplanet_df.columns:
            exoplanet_df['ra'] = exoplanet_df['sky_coord.ra']
            exoplanet_df['dec'] = exoplanet_df['sky_coord.dec']
        else:
            logging.warning(f"RA/DEC columns 'sky_coord.ra' and 'sky_coord.dec' not found. Exoplanet data might be incomplete.")
            return None

        exoplanet_df.to_csv(exoplanet_file, index=False)
        logging.info(f"Exoplanet data saved to {exoplanet_file}")

    logging.info(f"Number of exoplanets saved in file {exoplanet_filename}: {len(exoplanet_df)}")
    return exoplanet_df


####################################################################################
# download_star_data
# Download Star/Stellar data from TIC (TESS Input Catalog) catalog
####################################################################################

# Download star data for the specified region and save locally
# def download_star_data(ra_center, dec_center, search_radius):
# function definition allows to call this function without passing values for any/all of these parameters
def download_star_data(ra_center=None, dec_center=None, search_radius=None, tic_id=None):
    
    # Determine file naming based on whether a TIC ID or RA/Dec is used
    if tic_id:
        star_filename = f"star_tic_{tic_id}.csv"
    else:
        star_filename = f"star_ra{ra_center}_dec{dec_center}_radius{search_radius}.csv"

    star_file = os.path.join(BASE_DATA_DIR, star_filename)

    if os.path.exists(star_file):
        logging.info(f"Star data for the region already exists locally.")
        star_df = pd.read_csv(star_file)
    else:
        # If TIC ID is provided, query by TIC ID directly
        if tic_id:
            logging.info(f"Querying TIC catalog for TIC ID: {tic_id}")
            try:
                tic_data = Catalogs.query_object(f"TIC {tic_id}", catalog="TIC")
            except Exception as e:
                logging.error(f"An error occurred during TIC Catalog query for TIC ID {tic_id}: {e}")
                return None
        # Otherwise, query by RA/Dec and radius
        else:
            if ra_center is None or dec_center is None or search_radius is None:
                logging.error("RA, Dec, and search radius must be provided if TIC ID is not specified.")
                return None

            coord = SkyCoord(ra=ra_center * u.deg, dec=dec_center * u.deg)
            radius = search_radius * u.deg
            logging.info(f"Querying TIC catalog for stars in the region.")
            try:
                tic_data = Catalogs.query_criteria(
                catalog='TIC',
                radius=radius,
                coordinates=coord,
                Tmag=[TMAG_MIN, TMAG_MAX],  # Filter by TESS Magnitude range
                # plx=[1 / MAX_DISTANCE_PC, None],  # Filter by parallax for distance (<1000 pc)
                rad=[RAD_MIN, RAD_MAX],  # Filter by stellar radius
                Teff=[TEFF_MIN, TEFF_MAX],  # Filter by effective temperature
                # lumclass='DWARF',  # Filter for main-sequence stars (optional)
                objType='STAR'  # Only query stars
            )
            except Exception as e:
                logging.error(f"An error occurred during TIC Catalog query: {e}")
                return None

        star_df = tic_data.to_pandas()
        logging.info(f"Number of stars found: {len(star_df)}") 

        # Storing filtered star data to CSV    
        star_df.to_csv(star_file, index=False)
        logging.info(f"Star data saved to {star_file}")

    logging.info(f"Number of stars saved in file {star_filename}: {len(star_df)}")
    logging.info(f"Listing all the columns: {(star_df.columns)}")
    return star_df


############################################################################
# prepare_dataset
# Prepare the final dataset by combining exoplanet and star data
# This custom dataset will be used in the UI
############################################################################

def prepare_dataset(ra_center, dec_center, search_radius):
    final_dataset_filename = f"final_dataset_ra{ra_center}_dec{dec_center}_radius{search_radius}.csv"
    final_data_file = os.path.join(BASE_DATA_DIR, final_dataset_filename)

    # Check if the final dataset file exists
    if os.path.exists(final_data_file):
        logging.info(f"Final dataset for the region already exists.")
        star_df = pd.read_csv(final_data_file)
    else:
        # Proceed with dataset creation if it doesn't exist
        exoplanet_df = download_exoplanet_data(ra_center, dec_center, search_radius)
        star_df = download_star_data(ra_center, dec_center, search_radius)

        # Handle cases where data might be missing
        if exoplanet_df is None or star_df is None or exoplanet_df.empty or star_df.empty:
            logging.warning(f"Aborting final dataset creation - No valid data for RA={ra_center}, Dec={dec_center}.")
            return None

        # Ensure 'ra' and 'dec' columns are present in both dataframes
        if 'ra' not in exoplanet_df.columns or 'dec' not in exoplanet_df.columns:
            logging.error("Exoplanet data is missing 'ra' or 'dec' columns.")
            return None
        if 'ra' not in star_df.columns or 'dec' not in star_df.columns:
            logging.error("Star data is missing 'ra' or 'dec' columns.")
            return None

        # Call the cross-match method
        idx, sep_constraint, matched_exoplanets = cross_match_coordinates(exoplanet_df, star_df)

        # Assign c_has_exoplanet flag and add exoplanet info - custom calculated columns with 'c_' prefix
        star_df['c_has_exoplanet'] = sep_constraint 
        star_df['c_exoplanet_name'] = None 
        star_df['c_pl_orbper'] = None   # time it takes to complete one orbit around its star (data from exoplanet archive)
        star_df['c_pl_rade'] = None  # exoplanet radius (data from exoplanet archive)

        # Update only rows where 'c_has_exoplanet' is True 
        matched_exoplanets = exoplanet_df.iloc[idx[sep_constraint]]
        star_df.loc[sep_constraint, 'c_exoplanet_name'] = matched_exoplanets['pl_name'].values
        star_df.loc[sep_constraint, 'c_pl_orbper'] = matched_exoplanets['pl_orbper'].values
        star_df.loc[sep_constraint, 'c_pl_rade'] = matched_exoplanets['pl_rade'].values

        # Assign a spestar type classification based on effective temperature (Teff) ranges (O, B, A, F, G, K, M)
        logging.info(f"Setting up 'c_star_type' for all the stars as TIC catalog doesn't give star type (O,B,A,F,G,K,K)")  
        star_df['c_star_type'] = star_df['Teff'].apply(classify_star)

        # Save the final dataset
        star_df.to_csv(final_data_file, index=False)
        logging.info(f"Final star dataset saved to {final_data_file}")

        # Print count of stars with exoplanets and relevant exoplanet information
        # exoplanet_hosts = star_df[star_df['c_has_exoplanet'] == True]
        # print(f"Count of stars with exoplanets: {len(exoplanet_hosts)}")
        # for index, row in exoplanet_hosts.iterrows():
        #     logging.info(
        #         f"Host Star TIC ID (ID)        : {row['ID']}\n"
        #         f"Exoplanet Name (c_exoplanet_name)  : {row['c_exoplanet_name']}\n"
        #         f"Orbital Period (c_pl_orbper)   : {row['c_pl_orbper']} days\n"
        #         f"Planet Radius (c_pl_rade)      : {row['c_pl_rade']} Earth radii\n"
        #         f"RA (ra)                      : {row['ra']}\n"
        #         f"DEC (dec)                    : {row['dec']}\n"
        #         f"Tmag (Tmag)                  : {row['Tmag']}\n"
        #         f"Parallax (plx)               : {row['plx']}\n"
        #         f"Star Radius (rad)            : {row['rad']} solar radii\n"
        #         f"Temperature (Teff)           : {row['Teff']} K\n"
        #         f"Luminosity Class (lumclass)  : {row['lumclass']}\n"
        #         f"{'-' * 50}"
        #     )           

    exoplanet_hosts = star_df[star_df['c_has_exoplanet'] == True]
    logging.info(f"Count of stars with exoplanets in final star dataset: {len(exoplanet_hosts)}")
    logging.info(f"Final star dataset exists with row count= {len(star_df)}")     
    return star_df


######################################
# analyze_light_curve
# Function for Light Curve Analysis
######################################

def analyze_light_curve(tic_id):

    # Directory to store light curve data
    LIGHTCURVE_DATA_DIR = os.path.join(BASE_DATA_DIR, 'lightcurve_data')
    os.makedirs(LIGHTCURVE_DATA_DIR, exist_ok=True)

    # Option to limit sectors (e.g., sectors with good data)
    # sectors_to_use = [1, 2, 3, 4, 5, 6, 7, 8, 9]  # Adjust as needed
    sectors_to_use = range(1, 13)  # Adjust as needed

    # Loop through and return the result for the first successful sector (period, depth, planet_radius, and the phase/flux values)
    for sector in sectors_to_use:
        logging.info(f"Sector#{sector}")

        # Filename for the light curve data
        filename = f"lightcurve_TIC{tic_id}_sector{sector}.fits"
        filepath = os.path.join(LIGHTCURVE_DATA_DIR, filename)

        # Check if the light curve file exists locally
        if os.path.exists(filepath):
            logging.info(f"Loading light curve from {filepath}")
            lc = lk.read(filepath)

            # Comment this check later - After loading, the 'flux' attribute should be present 
            # so no need to check for 'pdcsap_flux' or 'sap_flux' columns
            if hasattr(lc, 'flux') and hasattr(lc, 'flux_err'):
                logging.info("Loaded light curve has 'flux' and 'flux_err' attributes.")
            else:
                logging.error("Loaded light curve does not have 'flux' or 'flux_err' attributes.")
                continue  # Try next sector
        else:
            logging.info(f"Downloading light curve for TIC {tic_id}, sector {sector}")
            # Search for SPOC light curve files with 2-minute cadence    
            search_result = lk.search_lightcurve(
                f"TIC {tic_id}",
                mission='TESS',
                author='SPOC',
                exptime=120,
                sector=sector
            )

            if len(search_result) == 0:
                logging.error(f"No SPOC light curve files found for TIC ID {tic_id}")
                continue  # Try next sector
            else:
                logging.info(f"Total search results = {len(search_result)}")
                logging.info(search_result)

            # Download and stitch the light curves
            try:
                logging.info(f"....Downloading and stitching light curves")
                lc_collection = search_result.download_all()
                lc = lc_collection.stitch()

                # Check available columns
                available_flux_columns = lc.columns
                logging.info(f"Available columns in light curve: {available_flux_columns}")

                # Prioritize pdcsap_flux, fallback to sap_flux
                if 'pdcsap_flux' in available_flux_columns:
                    flux = lc['pdcsap_flux']
                    flux_err = lc['pdcsap_flux_err']
                elif 'sap_flux' in available_flux_columns:
                    flux = lc['sap_flux']
                    flux_err = lc['sap_flux_err']
                    logging.warning("pdcsap_flux not available. Using sap_flux instead.")
                else:
                    logging.error("No valid flux column ('pdcsap_flux' or 'sap_flux') available.")
                    continue  # Try next sector

                # Update the light curve with the selected flux and flux_err
                lc = lc.copy()
                lc.flux = flux
                lc.flux_err = flux_err
        
                # Save the stitched light curve as a FITS file 
                lc.to_fits(filepath, overwrite=True) 
                logging.info(f"Saved light curve to {filepath}")

            except Exception as e:
                logging.error(f"Error downloading or stitching light curves: {e}")
                continue  # Try next sector

            # Before saving the file, verify that it contains the required columns -'pdcsap_flux' or 'sap_flux'
            if 'pdcsap_flux' in lc.columns or 'sap_flux' in lc.columns:
                # Proceed with saving
                logging.info("Required flux columns are available. Proceeding with save.")
            else:
                logging.error("Required flux columns ('pdcsap_flux' or 'sap_flux') are missing. Aborting save.")
                continue  # Try next sector  # Skip saving or handle the error as needed

        
        # At this point, we have lc (either from reading local file or querying lk.search_lightcurve) 
        # At this point lc is guaranteed to have 'flux' and 'flux_err' attributes needed for Plotly plotting 

        # Check available columns
        logging.info(f"Available columns in light curve: {lc.columns}")
        
        # Convert flux and flux_err to unitless numpy arrays
        lc.flux = lc.flux.value  # Remove units from flux
        lc.flux_err = lc.flux_err.value  # Remove units from flux_err
        
        ### Process Valid Light Curve Data:
        ##  If valid flux data is found, the code 
        ##      processes the light curve, removes outliers, flattens the curve, and performs a Box Least Squares (BLS) analysis to detect transits.
        ##  After analyzing the data, return the analysis results for the first successful sector and exits the loop.

        # Clean the light curve
        logging.info(f"....Cleaning light curve")
        lc = lc.remove_nans().remove_outliers()

        # Flatten the light curve to remove trends
        logging.info(f"....Flattening light curve")
        flat_lc = lc.flatten()

        # Use the Box Least Squares (BLS) method to search for transits
        logging.info(f"....Using BLS method to search for transits")     
        # Define the period grid explicitly
        periods = np.linspace(0.5, 10, 10000)  # Periods from 0.5 to 10 days
        periodogram = flat_lc.to_periodogram(
            method='bls',
            period=periods  # define a period grid from 0.5 to 10 days
        )

        best_period = periodogram.period_at_max_power
        transit_time = periodogram.transit_time_at_max_power
        depth = periodogram.depth_at_max_power

        # Calculate the planet radius
        star_info = Catalogs.query_object(f"TIC {tic_id}", catalog="TIC")[0]
        star_radius = star_info['rad']  # Stellar radius in solar radii
        if star_radius is None:
            logging.error(f"Stellar radius not available for TIC ID {tic_id}")
            continue  # Try next sector

        transit_depth = depth.value
        planet_radius = star_radius * np.sqrt(transit_depth)  # Planet radius in solar radii

        # Prepare data for plotting
        folded_lc = flat_lc.fold(period=best_period, epoch_time=transit_time)

        json_lightcurve_info = {
            'period': best_period.value,
            'depth': depth.value,
            'planet_radius': planet_radius,
            'sector': flat_lc.meta.get('SECTOR'),

            ## Updated plotting using Plotly: Folded light curve and time series light curve             
            # plotting 'folded light curve' needs folded_lc phase and flux data
            'phase': folded_lc.phase.value.tolist(),  
            'flux_phase': folded_lc.flux.value.tolist(),

            # plotting 'time series light curve' needs flat_lc time and flux data
            'time': flat_lc.time.value.tolist(),
            'flux_time': flat_lc.flux.value.tolist()
        }

        # Return Combined JSON response - light curve data and info box data
        return {
            'meta_observation_details': meta_observation_details(lc),
            'meta_stellar_details': meta_stellar_details(tic_id),
            'light_curve': json_lightcurve_info
        }
    
        
    # If all sectors fail
    logging.error(f"Analysis failed for TIC ID {tic_id}. No valid data in the specified sectors.")
    return None


######################################
# Main method
# Execution starts from here 
######################################

if __name__ == '__main__':
    logging.info(f"Processing dataset")

    ## #1 - Uncomment the following line to test datasets download (exoplanet, star and final dataset)
    # prepare_dataset(RA_CENTER, DEC_CENTER, SEARCH_RADIUS)  
    # download_star_data(tic_id=25155310)

    ## #2 - Uncomment the following lines to test folded light curve data using matplotlib (i.e. .png file will be created and downloaded)
    ## Below test will not work as matplotlib implementaiton is commented in all the functions above. 
    # logging.info("Testing analyze_light_curve function locally.")
    # test_tic_id = 25155310  # Replace with your TIC ID
    # result = analyze_light_curve(test_tic_id)
    # if result:
    #     print(f"Orbital Period: {result['period']:.5f} days")
    #     print(f"Transit Depth: {result['depth'] * 100:.6f}%")
    #     print(f"Estimated Planet Radius: {result['planet_radius']:.5f} Solar Radii")
    #     # Save the plot to a file for verification
    #     with open('light_curve_plot.png', 'wb') as f:
    #         f.write(base64.b64decode(result['plot_data']))
    #     print("Light curve plot saved as 'light_curve_plot.png'.")
    # else:
    #     print("Light curve analysis failed.")


    ## #3 - Uncomment the following lines to test folded light curve data using matplotlib (i.e. .png file will be created and downloaded)
    #### TBD - write test for folded light curve and time-series light curve function


    ## #4 - Uncomment Uncomment the following line when testing with frontend  
    app.run(debug=True)
