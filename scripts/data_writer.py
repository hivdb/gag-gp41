import os
import csv


def csv_writer(path, data, headers=None, writer_options=None):
    dirname = os.path.dirname(path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    if headers is None:
        headers, data = data
    writer_options = writer_options or {}
    with open(path, 'w') as fp:
        writer = csv.DictWriter(fp, headers, **writer_options)
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
