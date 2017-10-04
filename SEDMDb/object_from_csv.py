from __future__ import print_function
import csv
import SedmDb
from SedmDb_tools import dec_to_decimal, ra_to_decimal
import argparse

db = SedmDb.SedmDB(dbname='sedmdbtest', host='localhost')


def objects_from_csv(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile, )
        next(reader, None)
        for row in reader:  # For each row, create the dictionary and submit
            pardic = {'name': row[0], 'ra': ra_to_decimal(row[1]), 'dec': dec_to_decimal(row[2]), 'mag': row[3]}
            pardic['typedesig'] = 'f'
            db.add_object(pardic)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str)
    args = parser.parse_args()
    f = args.file
    objects_from_csv(f)
