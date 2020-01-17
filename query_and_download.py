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
from tqdm import tqdm
from joblib import Parallel, delayed
from astroquery.simbad import Simbad

from .query_object import query_object as _query_object

# --- setup a custom simbad query with the relevant fields
cSimbad = Simbad()
cSimbad.remove_votable_fields('coordinates')
cSimbad.add_votable_fields('otype(3)', 'z_value', 'ra(d)', 'dec(d)')


def query_and_download(database=None, querytype="both", inname=None,
                       inra=None, indec=None, verbose=False,
                       sdss_dr=15, sdss_onlyprim=True, conrad=3, outfolder=None,
                       parallel=1, overwrite=False, istart=0):
    """
    Main routine to query an online database through astroquery for a provided
    list of coordinates or names in a parellelized way and write the output
    results into a folder with one file per object. The feed_results_into_table
    routine does then gather the results into a table

    Parameters
    ----------
    database : TYPE, optional
        DESCRIPTION. The default is None. Name string of the database to query
        (SIMBAD, NED, SDSS)
    querytype : TYPE, optional
        DESCRIPTION. The default is "both". Type of query: by 'name', or
        'coordinates' or 'both' (first by name and if it does not return any
        results then by coordinates)
    inname : TYPE, optional
        DESCRIPTION. The default is None. Input object name list to query for.
    inra : TYPE, optional
        DESCRIPTION. The default is None. Input object coordinar RA to query
        for.
    indec : TYPE, optional
        DESCRIPTION. The default is None.  Input object coordinar DEC to query
        for.
    sdss_dr : TYPE, optional.
        DESCRIPTION. The default is 15. Data release number of SDSS to query
    sdss_onlyprim : TYPE, optional
        DESCRIPTION. The default is True. Flag to select whether to query only
        for primary science objects in SDSS or for all (including dublicates)
    conrad : TYPE, optional
        DESCRIPTION. The default is 3. Cone radius for the crossmatching by
        coordinates in arcsec
    outfolder : TYPE, optional
        DESCRIPTION. The default is None. Output folder
    parallel : TYPE, optional
        DESCRIPTION. The default is 1. Level of parallelisation
    overwrite : TYPE, optional
        DESCRIPTION. The default is False. Flag whether to overwrite previous
        output
    istart : TYPE, optional
        DESCRIPTION. The default is 0. Minimum starting index in the object
        list for the query (in case one does not want to restart at 0)

    Returns
    -------
    None.

    """


#    nname = 0
#    ncoord = 0
#    n_not_found = 0

    if inname is not None:
        n = len(inname)
    else:
        n = len(inra)
#
#    outsep = np.ma.zeros(n, dtype=float)
#    outname = np.ma.zeros(n, dtype=object)
#    outra = np.ma.zeros(n, dtype=float)
#    outdec = np.ma.zeros(n, dtype=float)
#    outz = np.ma.zeros(n, dtype=float)
#    outtype = np.ma.zeros(n, dtype=object)


    if database == "SDSS":

        # --- only the sciencePrimary?
        if sdss_onlyprim:
            sdsstable = "SpecObj"
        else:
            sdsstable = "SpecObjAll"

    else:
        sdsstable = None


    if verbose:
        print("Search parameters:")
        print("   - cone radius: ", conrad)
        print("   - query type: ", querytype)
        print("   - database: ", database)
        if database == "SIMBAD":
            print("   - selected fields: ", cSimbad.get_votable_fields())


    outfile = np.zeros(n, dtype=object)
    for i in range(n):
        outfile[i] = (outfolder + "/" + '{:.6f}'.format(inra[i])
                      + "_"+ '{:.6f}'.format(indec[i]) + "_"
                      + inname[i].replace(" ","") + ".csv")



    # --- loop over all the object
    if parallel == 1:
#        for i in tqdm(range(n)):
        for i in tqdm(np.arange(istart,n)):

            if verbose:
                print(i, inname[i], inra[i], indec[i])

            _ = _query_object(outfile[i], database=database, querytype=querytype,
                               inname=inname[i], inra=inra[i], indec=indec[i],
                               verbose=verbose, sdss_dr=sdss_dr,
                               sdss_onlyprim=sdss_onlyprim, conrad=conrad,
                               sdsstable=sdsstable, overwrite=overwrite)

#            return(res)

    else:
#    f = open(outfile, "w")

        _ = Parallel(n_jobs=parallel)(delayed(_query_object)(outfile[i], database=database, querytype=querytype,
                               inname=inname[i], inra=inra[i], indec=indec[i],
                               verbose=verbose, sdss_dr=sdss_dr,
                               sdss_onlyprim=sdss_onlyprim, conrad=conrad,
                               sdsstable=sdsstable, overwrite=overwrite) for i in np.arange(istart,n))

