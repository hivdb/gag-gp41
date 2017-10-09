import os
import csv


def csv_writer(path, data, headers=None):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    if headers is None:
        headers, data = data
    with open(path, 'w') as fp:
        writer = csv.DictWriter(fp, headers)
        writer.writeheader()
        writer.writerows(data)
    print('{} created'.format(path))


def data_writer(path, data):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(path, 'w') as fp:
        fp.write(data)
    print('{} created'.format(path))
