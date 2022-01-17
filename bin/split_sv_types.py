#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen


class VcfWriteProc:
    def __init__(self, pref, header, as_bcf):
        self.as_bcf = as_bcf
        self.fn = (pref + '.bcf') if as_bcf else (pref + '.vcf.gz')
        mode = '-Ob' if as_bcf else '-Oz'
        self.Popen = Popen(['bcftools', 'view', mode, '--no-version', '-o', self.fn],
                           stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.Popen.stdin, 'wu', header=header)

    def index(self):
        flag = '' if self.as_bcf else '-t'
        popen = Popen(['bcftools', 'index', flag, self.fn],
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


def get_sv_type(record):
    if 'SVTYPE' in record.info:
        return record.info.get('SVTYPE')
    return 'NONE'


def main(input, pref, as_bcf, do_index):
    vf_in = VariantFile(input)
    vf_outs = {}
    for rec in vf_in:
        sv_type = get_sv_type(rec)
        if sv_type not in vf_outs:
            vf_outs[sv_type] = VcfWriteProc("{}.{}".format(pref, sv_type),
                                            vf_in.header, as_bcf)
        vf_outs[sv_type].write(rec)
    for vf in vf_outs.values():
        vf.close(index=do_index)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='input vcf/vcf.gz/bcf file')
    parser.add_argument('--pref', default='output', help='output variant files prefix')
    parser.add_argument('--bcf', action='store_true', help='write bcf instead of vcf.gz')
    parser.add_argument('--index', action='store_true', help='index output files')
    args = parser.parse_args()
    main(args.input, args.pref, args.bcf, args.index)
