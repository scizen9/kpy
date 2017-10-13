from __future__ import print_function
import csv
import SedmDb
from SedmDb_tools import dec_to_decimal, ra_to_decimal
import argparse

db = SedmDb.SedmDB(dbname='sedmdb', host='localhost')


def objects_from_csv(filename):
    """
    Opens a csv file with format name, ra, dec, mag(optional) and inputs it into the database
    NOTE: assumes fixed object
    """
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader, None)
        for row in reader:  # For each row, create the dictionary and submit
            pardic = {'name': row[0], 'ra': ra_to_decimal(row[1]), 'dec': dec_to_decimal(row[2])}
            pardic['typedesig'] = 'f'
            if len(row) > 3:
                mag = row[3] #not currently stored in object table
            obj = db.add_object(pardic)
            print(obj)


def objects_from_csv_allow_dup(filename):
    with open(filename, 'r') as csvfile:
        reader = csv.reader(csvfile, )
        next(reader, None)
        for row in reader:  # For each row, create the dictionary and submit
            db.execute_sql("INSERT INTO object (name, ra, dec, typedesig) VALUES ('%s', %s, %s, '%s')"
                           % (row[0], ra_to_decimal(row[1]), dec_to_decimal(row[2]), 'f'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str)
    args = parser.parse_args()
    if args.file:
        f = args.file
        objects_from_csv_allow_dup(f)
