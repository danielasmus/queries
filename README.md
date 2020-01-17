Collection of python routines to query astronomical databases

The main routine is query_and_download which calls query_object.
feed_results_into_table is then called after to feed results of the query_and_download into one table.

Used for data collection of LASr I (Asmus et al. 2020)

Written in Python 3 (might not run in Python 2)

Requires: numpy, joblib, tqdm, scipy, astropy, astroquery.

Created by Daniel Asmus


