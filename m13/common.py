import matplotlib.pyplot as plt
import seaborn as sns
import astropy.units as u
import numpy as np


def m13_info():
    return { 'ra':250.423475*u.deg, 
            'dec':36.46131944*u.deg, 
            'pmra':-3.18*u.mas/u.year, 
            'pmdec': -2.56*u.mas/u.year,
            'diameter': 33.0*u.arcmin,
            'parallax': 0.127*u.mas, 
            'age': 11.55*u.Myr,
            'distance': 6.8*u.kpc,
            'reddening': 0.04,
    }

import read_mist_models
from astropy.table import Table
from astropy.coordinates import Distance
def load_isochrone(filename):
    
    """ Read the MIST file into an isochrone table 
        Only includes the main sequence """
    iso = read_mist_models.ISOCMD(filename)
    iso_array = iso.isocmds[0]

    # select stars in main sequence and red giant phases 0<=phase<3
    phase = iso_array['phase']
    phase_mask = (phase>=0) & (phase<3)
    iso_main_sequence = iso_array[phase_mask]
    mag_g = iso_main_sequence['PS_g']
    color_gi = iso_main_sequence['PS_g'] - iso_main_sequence['PS_i']

    isochrone = Table([mag_g, color_gi], names=('mag_g', 'color_gi'))
    return isochrone

def adjust_isochrone(isochrone, distance, reddening, extinction):
    """ Adjust the color and magnitude values of this isochrone based on the given parameters
      Returns a new deep copy of the isochrone table  """
    
    # make a deep copy of the isochrone table
    isochrone = isochrone.copy(copy_data=True)

    # get a distance modulus from the distance value
    distance_obj = Distance(distance)
    distmod =  distance_obj.distmod.value

    # adjust the isochrone values based on distance, reddening, etc 
    isochrone['mag_g'] += distmod
    isochrone['color_gi'] += reddening

    #need to take into account extinction
    #mag_g += 0.12
    return isochrone


def polygon (a, b):
    return np.append(a, b[::-1])

def plot_cmd(stars):
    y=stars['g_mean_psf_mag']
    x=stars['g_mean_psf_mag'] - stars['i_mean_psf_mag']

    plt.plot(x,y,'ko', markersize=.1, alpha=.5)
    
    plt.xlim(0,3)
    plt.ylim(14,22)
    plt.gca().invert_yaxis()
    
    plt.ylabel('Magnitude $(g)$')
    plt.xlabel('Color $(g-i)$')

def plot_number_density(stars):
    x=stars['ra']
    y=stars['dec']
    f, ax = plt.subplots(figsize=(6, 6))
    #sns.histplot(x=x, y=y, bins=100, pthresh=.1)
    sns.kdeplot(x=x, y=y, levels=40, linewidths=1)

    plt.xlim(249, 252)
    plt.ylim(35, 38)

def plot_scatter(stars):
    x=stars['ra']
    y=stars['dec']
    f, ax = plt.subplots(figsize=(6, 6))
    #sns.histplot(x=x, y=y, bins=100, pthresh=.1)
    plt.plot(x, y, 'ko', markersize=.1, alpha=1)

    plt.xlim(249, 252)
    plt.ylim(35, 38)




def create_isochrone_polygon(isochrone, mag_range, color_range):
    """ Create a polygon to select a region around the isochrone in a color-magnitude diagram """

    # break out the range parameters
    (mag_top, mag_bottom) = mag_range
    (color_left, color_right) = color_range
    
    # filter to the magnitude range
    mag_g = isochrone['mag_g']
    g_mask = (mag_g > mag_top) & (mag_g <mag_bottom)
    iso_masked = isochrone[g_mask]
    mag_g = iso_masked['mag_g']

    # create a polygon to select a region around the isochrone
    color_gi = iso_masked['color_gi']
    color_lhs = color_gi + color_left
    color_rhs = color_gi + color_right

    # create a table with the polygon
    color_polygon = polygon(color_lhs, color_rhs)
    mag_polygon = polygon(mag_g, mag_g)
    return Table([color_polygon, mag_polygon], names=['color_gi', 'mag_g'])

from matplotlib.patches import Polygon
import pandas as pd

def filter_to_isochrone_polygon(stars, polygon):
    """ Filter the stars to the given color-magnitude polygon """

    # working in pandas
    poly_tester = Polygon(polygon.to_pandas())
    points = pd.DataFrame()

    # set up the dataframe of points 
    # ! note the polygon tester is sensitive to the order that these are added!
    points['color'] = stars['g_mean_psf_mag'] - stars['i_mean_psf_mag']
    points['mag'] = stars['g_mean_psf_mag']

    # mask to stars that are in the polygon
    inside = poly_tester.contains_points(points)
    iso_stars = stars[inside]
    return iso_stars