#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__version__ = "1.0.0"

"""
HISTORY:
    - 2020-01-17: created by Daniel Asmus


NOTES:
    -

TO-DO:
    - clean parallelisation
"""


import numpy as np
import os
import time
import datetime
from tqdm import tqdm
from astropy.io import ascii
import astropy.table as T
from astropy import units as u
from astropy import coordinates
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad


# --- setup a custom simbad query with the relevant fields
cSimbad = Simbad()
cSimbad.remove_votable_fields('coordinates')
cSimbad.add_votable_fields('otype(3)', 'z_value', 'ra(d)', 'dec(d)')


#%%

# --- helper routing
def timestamp():
    """
    return a handy string with the current time
    """

    return(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))


#%%
# --- helper routine
def ang_dist(r1, d1, r2, d2):
    """
    return the angular distance between two positions in the sky in arcsec
    positions cn be provided in degrees (float) or in sexagesimal units (string)
    """

    if np.isreal(r1):
        r1u = u.deg
    else:
        r1u = u.hourangle

    if np.isreal(r2):
        r2u = u.deg
    else:
        r2u = u.hourangle

    d1u = u.deg
    d2u = u.deg

    c1 = SkyCoord(r1, d1, unit=(r1u, d1u))
    c2 = SkyCoord(r2, d2, unit=(r2u, d2u))

#    print(r1u, d1u, r2u, d2u)
#    print(c1)
#    print(c2)

    d = c2.separation(c1)

    return(d.arcsec)


#%%
# --- helper routine
def read_res_file(d, database, ids, i, infile, inname, inra, indec, fail,
                  found_by_name=None, found_by_coord=None,
                  racol=None, deccol=None, typecol=None,
                  preftypes=None, verbose=False):

    """
    Helper routine to load the results of a query from a file and add them to a
    table
    """

    # --- check whether there is a matching file with the query results
    if not os.path.isfile(infile):

        fail[i] = 1

        if verbose:
            print(" - WARNING: No results found for object: ", i, inname)

        return(-1)


    else:

        res =  ascii.read(infile, header_start=0, delimiter=',', guess=False)

        nres = len(res)

        if nres == 0:
            fail[i] = 1

            if verbose:
                print(" - WARNING: Length of results = 0: ", i, inname)

            return(-1)

        # --- convert types into correct strings
        for j in range(nres):
            res[typecol][j] = res[typecol][j].replace("b'", "").replace("'", "")

        # --- the SDSS query at the moment is still not astroquery
        #     conform and does not return the separation so that
        #     we have to calculate it

        if "Separation" not in res.colnames:
            res.add_column(T.MaskedColumn(np.ma.zeros(nres, dtype=float),
                                          name="Separation"))

            if database != "SDSS":
                found_by_name[i] = 1
            else:
                found_by_coord[i] = 1

        # --- it seems to be added but empty in fact
        else:

            # --- if it exists, make sure the column has right type
            res['Separation'] =res['Separation'].astype(float)

            if  T.MaskedColumn(res['Separation']).mask[0]:
                found_by_name[i] = 1

            else:
                found_by_coord[i] = 1

        # --- just to be sure, (re)compute the separations
        for j in range(nres):
            res['Separation'][j] = ang_dist(inra, indec,
                                              res[racol][j],
                                              res[deccol][j])


        # --- sort the results by separation so that the closest comes first
        res = res[np.argsort(res['Separation'])]


        # --- get the objects in the result table with preferred type
        idg =  [x for x, t in enumerate(res[typecol])
                if t in preftypes]

        if len(idg) > 0:
            sel = idg[0]
        else:
            # --- if no preferred object is available, just take the closest
            sel = 0

            if verbose:
                print(" - WARNING: No object of preferred type found: ", i, inname[i])


        if verbose:
            print(inname, " ---> " ,res[sel])


        # --- fill the corresponding table values
        if database == "NED":

            d["NED_name"][ids[i]] = res['Object Name'][sel].replace("b'", "").replace("'", "")

            if not T.MaskedColumn(res["Redshift"]).mask[sel]:
                d["NED_redshift"][ids[i]] = float(res["Redshift"][sel])

            d["NED_RA_deg"][ids[i]] = float(res[racol][sel])
            d["NED_DEC_deg"][ids[i]] = float(res[deccol][sel])
            d["NED_type"][ids[i]] = res[typecol][sel]


        elif database == "SIMBAD":

            d["CDS_name"][ids[i]] = res['MAIN_ID'][sel]

            if "b'" in d["CDS_name"][ids[i]]:
                d["CDS_name"][ids[i]] = d["CDS_name"][ids[i]].replace("b'", "").replace("'", "")
            elif 'b"' in d["CDS_name"][ids[i]]:
                d["CDS_name"][ids[i]] = d["CDS_name"][ids[i]].replace('"b""', "").replace('"""', "")

            if not T.MaskedColumn(res["Z_VALUE"]).mask[sel]:
                d["CDS_redshift"][ids[i]] = float(res["Z_VALUE"][sel])

            d["CDS_RA_deg"][ids[i]] = float(res[racol][sel])
            d["CDS_DEC_deg"][ids[i]] = float(res[deccol][sel])
            d["CDS_type"][ids[i]]  = res[typecol][sel]


        elif database == "SDSS":

            # --- construct the SDSS name
            coords = coordinates.SkyCoord(ra=res[racol][sel]*u.degree,
                                          dec=res[deccol][sel]*u.degree)

            d["SDSS_name"][ids[i]] = 'SDSS J{0}{1}'.format(coords.ra.to_string(
                                                            unit=u.hourangle,
                                                            sep='',
                                                            precision=2,
                                                            pad=True),
                                                           coords.dec.to_string(
                                                            sep='',
                                                            precision=1,
                                                            alwayssign=True,
                                                            pad=True))

            if not T.MaskedColumn(res["z"]).mask[sel]:
                d["SDSS_redshift"][ids[i]] = float(res["z"][sel])

            d["SDSS_RA_deg"][ids[i]] = float(res[racol][sel])
            d["SDSS_DEC_deg"][ids[i]] = float(res[deccol][sel])
            d["SDSS_spectype"][ids[i]] = res[typecol][sel]
            d["SDSS_specID"][ids[i]] = res['specobjid'][sel]
            if not T.MaskedColumn(res['subclass']).mask[sel]:
                d["SDSS_specclass"][ids[i]] = res['subclass'][sel].replace("b'", "").replace("'", "")
            d["SDSS_warn"][ids[i]] = res['zWarning'][sel]
            d["SDSS_redshift_unc"][ids[i]] = res['zErr'][sel]

            if res['type'][sel] == 3:
                d["SDSS_phottype"][ids[i]] = "GALAXY"

            elif res['type'][sel] == 6:
                d["SDSS_phottype"][ids[i]] = "STAR"


        # --- associate the right coordinate offsets
        if database == "NED":
            if d["origin"][ids[i]] == "SIMBAD":
                d["NED-CDS_sep_as"][ids[i]] = res['Separation'][sel]
            elif d["origin"][ids[i]] == "SDSS":
                d["NED-SDSS_sep_as"][ids[i]] = res['Separation'][sel]
            elif d["origin"][ids[i]] == "2MRS":
                d["NED-2MRS_sep_as"][ids[i]] = res['Separation'][sel]
            elif d["origin"][ids[i]] == "B70":
                d["NED-B70_sep_as"][ids[i]] = res['Separation'][sel]

        elif database == "SIMBAD":
            if d["origin"][ids[i]] == "NED":
                d["NED-CDS_sep_as"][ids[i]] = res['Separation'][sel]

        elif database == "SDSS":
            if d["origin"][ids[i]] == "NED":
                d["NED-SDSS_sep_as"][ids[i]] = res['Separation'][sel]


    return(0)



#%%
# --- main routine

def feed_results_into_table(d, ids, infolder, innames, inra, indec, database,
                            verbose=None, preftypes=None, parallel=1):

    """
    Take the results of the query routine in this package and compile them into
    a table. Parallelisation does not work here unfortunately.
    """

    nids = len(ids)

    nfail = 0
    n_found_by_name = 0
    n_found_by_coord = 0

    if database == "NED":

        typecol = 'Type'
        racol = 'RA(deg)'
        deccol = 'DEC(deg)'

        if preftypes is None:
            preftypes = ['G', 'QSO']


    elif database == "SIMBAD":

        # --- add the redshift as an output column
        cSimbad = Simbad()
        cSimbad.add_votable_fields('otype(3)', 'z_value', 'ra(d)', 'dec(d)')

        racol = 'RA_d'
        deccol =  'DEC_d'
        typecol = 'OTYPE_3'

        if preftypes is None:
            preftypes = ['AGN', 'LIN', 'SyG', 'Sy1', 'Sy2', 'Bla',
                                     'BLL', 'OVV', 'QSO', 'AG?', 'Q?', 'Bz?',
                                     'BL?', 'EmG', 'SBG', 'bCG', 'G', 'GiC',
                                     'BiC', 'GiG', 'GiP', 'rG', 'H2G', 'LSB',
                                     'IG', 'PaG', 'G?']


    elif database == "SDSS":

        racol = 'ra'
        deccol = 'dec'
        typecol = 'class'

        if preftypes is None:
            preftypes = ['GALAXY', 'QSO']


    if verbose:
        print("   - database: ", database)
        print("   - racol: ", racol)
        print("   - deccol: ", deccol)


    tstart = time.time()
    print(timestamp(), " -- Start")

    infiles = np.zeros(nids, dtype=object)

    fail = np.zeros(nids, dtype=int)
    found_by_name = np.zeros(nids, dtype=int)
    found_by_coord = np.zeros(nids, dtype=int)

    for i in range(nids):

#        infiles[i] = (infolder + "/Qno" + str(i) + "_"
#                      + innames[i].replace(" ","") + ".csv")

        infiles[i] = (infolder + "/" + '{:.6f}'.format(inra[i])
                      + "_"+ '{:.6f}'.format(indec[i]) + "_"
                      + innames[i].replace(" ","") + ".csv")


    # --- loop over all the searched objects
    if parallel == 1:
        for i in tqdm(range(nids)):
           _ =  read_res_file(d, database, ids, i, infiles[i], innames[i],
                              inra[i], indec[i], fail,
                              found_by_name=found_by_name,
                              found_by_coord=found_by_coord,
                              racol=racol, deccol=deccol, typecol=typecol,
                              preftypes=preftypes)

#    ---- parallel does not work currently
#    else:
#         _ = Parallel(n_jobs=parallel)(delayed(read_res_file)(d, database, ids,
#                              i, infiles[i], innames[i],
#                              inra[i], indec[i], fail,
#                              found_by_name=found_by_name,
#                              found_by_coord=found_by_coord,
#                              racol=racol, deccol=deccol, typecol=typecol,
#                              preftypes=preftypes) for i in range(nids))

    tend = time.time()
    print(timestamp(), " -- End; Duration: ", tend - tstart)

    n_found_by_name = np.sum(found_by_name)
    n_found_by_coord = np.sum(found_by_coord)
    nfail = np.sum(fail)


    print("Total number searched: ", nids)
    print(" - Found matches by name: ", n_found_by_name)
    print(" - Found matches by coord: ", n_found_by_coord)
    print(" - Objects not found or error: ", nfail)
    print(" - Check sum: ", n_found_by_name + n_found_by_coord + nfail)


    return(n_found_by_name, n_found_by_coord, nfail)

