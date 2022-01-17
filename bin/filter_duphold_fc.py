#!/usr/bin/env python3

from pysam import VariantFile
import argparse
import subprocess
from subprocess import Popen


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


def main(input, output, max_del_fc, min_dup_fc, min_len):
    vf_in = VariantFile(input)
    header = vf_in.header
    header.filters.add('BAD_FC', None, None, 'Average fold change too high for DEL or too low for DUP')
    vf_out = VcfWriteProc(output, header)
    for rec in vf_in:
        if not rec.filter.keys():
            rec.filter.add('PASS')
        sv_type = rec.info['SVTYPE']
        if (sv_type == 'DEL' or sv_type == 'DUP') and abs(rec.info['SVLEN'][0]) > min_len:
            non_ref = [any(s['GT']) for s in rec.samples.values()]
            fc = [s['DHFFC'] for s in rec.samples.values()]
            if any(non_ref):
                mean_fc = sum([fc for (fc, b) in zip(fc, non_ref) if b]) / sum(non_ref)
                keep = (mean_fc < max_del_fc) if sv_type == 'DEL' else (mean_fc > min_dup_fc)
                if not keep:
                    rec.filter.add('BAD_FC')
        vf_out.write(rec)

    vf_out.close(index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default='-', help='input variant file')
    parser.add_argument('--output', default='output', help='output variant filename')
    parser.add_argument('--max-del-fc', type=float, default=0.7, help='Max average FC for deletion call')
    parser.add_argument('--min-dup-fc', type=float, default=1.25, help='Min average FC for duplication')
    parser.add_argument('--min-len', type=int, default=500, help='Min SV length to apply FC filter')
    args = parser.parse_args()
    main(args.input, args.output, args.max_del_fc, args.min_dup_fc, args.min_len)
