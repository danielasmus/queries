#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.0"

"""
HISTORY:
    - 2020-01-15: created by Daniel Asmus


NOTES:
    -

TO-DO:
    -
"""


import numpy as np
import os
import time
from astropy.io import ascii
from astropy import units as u
from astropy import coordinates
from astroquery.ned import Ned
from astroquery.simbad import Simbad
from astroquery.sdss import SDSS

cSimbad = Simbad()
cSimbad.remove_votable_fields('coordinates')
cSimbad.add_votable_fields('otype(3)', 'z_value', 'ra(d)', 'dec(d)')


def query_object(outfile, database=None, querytype="both", inname=None,
                 inra=None, indec=None, verbose=False,
                 sdss_dr=15, sdss_onlyprim=True, conrad=3,
                 sdsstable=None, overwrite=False):

    """
    Query an online database for an object either by name or coordinates
    (or both) and write the results into an output file. Mostly this is a
    subroutine for the query_and_download routine which does the query for a
    list of objects/coordinates in a parallelized way. (called by query_and_download)
    """


    if not overwrite and os.path.isfile(outfile):
        if verbose:
            print("File exists! Skipping to next...")
        return(0)

    # --- first try to find a match by name
    if querytype == "both" or querytype == "name":

        if verbose:
            print("   - Searching by name ...")

        try:
            if database == "NED":
                res= Ned.query_object(inname)
            elif database == "SIMBAD":
                time.sleep(1)
                res = cSimbad.query_object(inname)

#            nname +=1
            _ = res.colnames

            if len(res) == 0:
                namefail = True
            else:
                namefail = False

            if verbose:
                if namefail:
                    print("   - object not found by name")
                else:
                    print("   - Found by name: ", res)

        except:
            namefail = True

            if verbose:
                print("   - object not found by name")

            if querytype == "name":
                if verbose:
                    print("Object not found: ", inname, inra,
                          indec)

                return(-1)

    # --- if coordinate search or by name fails, try by coordinates:
    if querytype == "coordinates" or (querytype == "both" and namefail):

        if verbose:
            print("   - Searching by coordinates ...")

        try:
            coord = coordinates.SkyCoord(ra=inra, dec=indec,
                                         unit=(u.deg, u.deg), frame='fk5')

            if verbose:
                print(coord)

            if database == "NED":
                res= Ned.query_region(coord, radius=conrad*u.arcsec,
                                      equinox='J2000')

                if len(res) == 0:
                    if verbose:
                        print("   - Object not found: ", inname, inra, indec)

                    return(-1)


            elif database == "SIMBAD":

                time.sleep(1)

                res = cSimbad.query_region(coord, radius=conrad*u.arcsec)

                print(" query done..")

            elif database == "SDSS":

                # --- query region only looks for coordinates in the photo
                #     catalog which for some reason is null for many
                #     galaxies. Therefore, we have to query by sql

                sqlquery = ("SELECT TOP 50 s.ra, s.dec, s.specobjid, s.z, \
                             s.zErr, s.zWarning, s.sciencePrimary, \
                             s.class, s.subclass, p.type FROM .."
                             + sdsstable + " as\
                             s LEFT JOIN ..PhotoObj AS p ON \
                             s.bestObjID=p.objID JOIN \
                             dbo.fGetNearby" + sdsstable + "Eq("
                             + str(inra) + "," + str(indec) + ","
                             + str(conrad/60.0) + ") \
                             AS b ON  b.SpecobjID = S.SpecobjID")

                res = SDSS.query_sql(sqlquery, data_release=sdss_dr)

#                    res = None
                if verbose:
                    print("   - Found by coordinates: ", res)


#                        return(res)

                    res = res[np.where(np.array(res['sciencePrimary']) == 1)[0]]

            if verbose:
                print("   - Found by coordinates: ", res)

        except:
            # --- if both fail, object is not in database apparently
    #            print("ERROR: No object found neither by name nor by coordinates!")

            if verbose:
                print("   - Object not found: ", inname, inra, indec)
            return(-1)

#    f.write(res)
    if res is not None:
        if len(res) > 0:
#        return(res)
            res.write(outfile, delimiter=',', format='ascii',
                      fill_values=[(ascii.masked, '')], overwrite=True)


    return(1)


