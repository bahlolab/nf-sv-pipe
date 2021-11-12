#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen
import gzip


class VcfWriteProc:
    def __init__(self, fn, header):
        self.fn = fn
        self.Popen = Popen(['bcftools', 'view', '-Oz', '--no-version', '-o', self.fn],
                           stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.Popen.stdin, 'wu', header=header)

    def index(self):
        popen = Popen(['bcftools', 'index', '-t', self.fn],
                      stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        return popen.poll()

    def close(self, index=False):
        self.variantFile.close()
        self.Popen.communicate()
        if index:
            self.Popen.poll()
            return self.index()
        return self.Popen.poll()

    def write(self, rec):
        self.variantFile.write(rec)


def read_ids(id_files):
    ids = {}
    for fn in id_files:
        with gzip.open(fn, 'rt') as handle:
            for id in handle:
                ids[id.strip()] = 1
    return ids


def main(out, id_files):
    ids = read_ids(id_files)
    vf_in = VariantFile('-')
    vf_out = VcfWriteProc(out, vf_in.header)
    for rec in vf_in:
        if any([x in ids for x in rec.info.get('IDLIST')]):
            rec.filter.clear()
            rec.filter.add('PASS')
            vf_out.write(rec)
    vf_out.close(index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='output', help='output variant filename')
    parser.add_argument('--ids', action='append', help='set of pass ids')
    args = parser.parse_args()
    main(args.out, args.ids)
