#!/usr/bin/env python3
import os.path

from pysam import VariantFile
import re
import argparse
import subprocess
from subprocess import Popen

MODE = {'bcf': 'b', 'vcf': 'v', 'gz': 'z'}
REMOVE = ['^INFO/SVTYPE', '^INFO/SVLEN', '^INFO/END', '^INFO/STRANDS',
          '^INFO/PRECISE', '^INFO/IMPRECISE', 'FORMAT']


class VcfWriter:
    def __init__(self, output, header, remove=None):
        mode = MODE[file_ext(output)]
        self.output = output
        if remove:
            cmd = ['bcftools', 'annotate', '--no-version', '--remove', ','.join(remove),
                   '-O{}'.format(mode), '-o', output]
        else:
            cmd = ['bcftools', 'view', '--no-version', '-O{}'.format(mode), '-o', output]
        self.proc = Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.variantFile = VariantFile(self.proc.stdin, 'wu', header=header)

    def index(self):
        proc = Popen(['bcftools', 'index', self.output],
                      stdin=subprocess.PIPE, stdout=subprocess.DEVNULL)
        return proc.poll()

    def close(self, index=False):
        self.variantFile.close()
        self.proc.communicate()
        if index:
            self.proc.poll()
            return self.index()
        return self.proc.poll()

    def write(self, rec):
        self.variantFile.write(rec)


def file_ext(filename):
    sp = filename.split('.')
    if len(sp) > 1:
        return sp[-1]


def get_sv_type(record):
    if 'SVTYPE' in record.info:
        return record.info.get('SVTYPE')
    return 'NONE'


def strip_vcfs(inputs, work, split_types):
    vcfs = {}
    for input in inputs:
        pref = '{}/{}'.format(work, re.sub('(.vcf(.gz)?)|(.bcf)$', '', input))
        vf = VariantFile(input)
        outs = {}
        for rec in vf:
            sv_type = get_sv_type(rec) if split_types else 'ALL'
            if sv_type not in outs:
                outs[sv_type] = VcfWriter("{}.{}.vcf".format(pref, sv_type), vf.header, remove=REMOVE)
            outs[sv_type].write(rec)
        for sv_type, vw in outs.items():
            if sv_type not in vcfs:
                vcfs[sv_type] = [vw.output]
            else:
                vcfs[sv_type].append(vw.output)
            vw.close()
    return vcfs


def jasmine(input_vcfs, out_vcf, jasmine_dir, work, threads=2, max_dist_linear=0.2, min_overlap=0.8):
    for sv_type in input_vcfs:
        input_list = '{}/{}.txt'.format(work, sv_type)
        output = '{}/{}.vcf'.format(work, sv_type)
        with open(input_list, "w") as out:
            out.write("\n".join(input_vcfs[sv_type]))
        cmd = ['java', '-cp', '{}/jasmine_iris.jar:{}/jasmine.jar'.format(jasmine_dir, jasmine_dir), 'Main',
               'iris_args=samtools_path=samtools,racon_path=racon,minimap_path=minimap2',
               'file_list={}'.format(input_list), 'out_file={}'.format(output), 'threads={}'.format(threads),
               'out_dir={}'.format(work), 'max_dist_linear={}'.format(max_dist_linear),
               'min_overlap={}'.format(min_overlap)]
        print(' '.join(cmd))
        proc = Popen(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.communicate()
        proc = Popen(['bcftools', 'view', output, '-Oz', '-o', '{}.gz'.format(output)],
                     stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.communicate()


def main(inputs, output, work, split_types, jasmine_jar):
    if not os.path.exists(work):
        os.makedirs(work)
    stripped = strip_vcfs(inputs, work, split_types)
    jasmine(stripped, output, jasmine_jar, work)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+', help='input sample vcfs to be merged')
    parser.add_argument('--output', default='output', help='output variant filename')
    parser.add_argument('--work', default='work', help='working directory')
    parser.add_argument('--jasmine-dir', default='/bin', help='path to jasmine jar file')
    parser.add_argument('--split-types', action='store_true', help='split sv types prior to jasmine process')
    args = parser.parse_args()
    main(args.inputs, args.output, args.work, args.split_types, args.jasmine_dir)
