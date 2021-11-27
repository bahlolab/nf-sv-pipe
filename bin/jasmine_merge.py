#!/usr/bin/env python3
import os.path
from pysam import VariantFile
import re
import argparse
import subprocess
from subprocess import Popen
import time

MODE = {'bcf': 'b', 'vcf': 'v', 'gz': 'z'}
REMOVE = ['^INFO/SVTYPE', '^INFO/SVLEN', '^INFO/END', '^INFO/STRANDS',
          '^INFO/PRECISE', '^INFO/IMPRECISE', 'FORMAT']
WORK = 'work'


class VcfWriter:
    def __init__(self, output, header, drop=False, remove=None, sort=False):
        mode = MODE[file_ext(output)]
        self.output = output
        cmd = []
        if remove:
            cmd.extend(['bcftools', 'annotate', '--no-version', '-Ou', '--remove', ','.join(remove), '|'])
        if sort:
            cmd.extend(['bcftools', 'sort', '-Ou', '-T', WORK, '|'])
        cmd.extend(['bcftools', 'view', '--no-version', '-O{}'.format(mode), '-o', output])
        if drop:
            cmd.extend(['-G'])
        self.proc = Popen(' '.join(cmd), stdin=subprocess.PIPE,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                          shell=True)
        self.variantFile = VariantFile(self.proc.stdin, 'wu', header=header)

    def index(self):
        proc = Popen(['bcftools', 'index', '-f', self.output],
                     stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.wait()
        return proc.poll()

    def close(self, index=False):
        self.variantFile.close()
        self.proc.communicate()
        if index:
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


def proc_vcfs(inputs, sv_types):
    vcfs = {}
    records = {}
    for input in inputs:
        pref = re.sub('(.vcf(.gz)?)|(.bcf)$', '', input)
        vf = VariantFile(input)
        # samples.extend(list(vf.header.samples))
        outs = {}
        records[pref] = {'HEADER': vf.header}
        for i, rec in enumerate(vf):
            sv_type = get_sv_type(rec)
            if sv_types and sv_type not in sv_types:
                continue
            rec.id = sv_type + '_' + pref + '_' + str(i)
            if sv_type not in outs:
                outs[sv_type] = VcfWriter("{}/{}.{}.vcf".format(WORK, pref, sv_type), vf.header, remove=REMOVE)
                records[pref][sv_type] = {}
            outs[sv_type].write(rec)
            records[pref][sv_type][rec.id] = rec
        for sv_type, vw in outs.items():
            if sv_type not in vcfs:
                vcfs[sv_type] = [vw.output]
            else:
                vcfs[sv_type].append(vw.output)
            vw.close()
    return vcfs, records


def get_jasmine_cp(jasmine_dir):
    jasmine_jar = jasmine_dir + '/jasmine.jar'
    iris_jar = jasmine_dir + '/jasmine_iris.jar'
    if not os.path.isfile(jasmine_jar) or not os.path.isfile(iris_jar):
        print("jasmine.jar/jasmine_iris.jar not found in {}".format(jasmine_dir))
        exit(1)
    return iris_jar + ':' + jasmine_jar


def jasmine(input_vcfs, jasmine_cp, args):
    vcfs = {}
    for sv_type in input_vcfs:
        input_list = '{}/{}.txt'.format(WORK, sv_type)
        output = '{}/{}.vcf'.format(WORK, sv_type)
        output_bcf = '{}/{}.bcf'.format(WORK, sv_type)
        with open(input_list, "w") as out:
            out.write("\n".join(input_vcfs[sv_type]))
        cmd = ['java', '-cp', jasmine_cp, 'Main',
               'iris_args=samtools_path=samtools,racon_path=racon,minimap_path=minimap2',
               '--keep_var_ids',
               'file_list={}'.format(input_list),
               'out_file={}'.format(output),
               'out_dir={}'.format(WORK)]
        cmd.extend(args)
        print(' '.join(cmd))
        proc = Popen(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.wait()
        if proc.poll() != 0:
            print("jasmine error")
            exit(1)
        sort_and_index(output, output_bcf)
        vcfs[sv_type] = output_bcf
    return vcfs


def sort_and_index(input, output):
    mode = MODE[file_ext(output)]
    proc = Popen(['bcftools', 'sort', input, '-O{}'.format(mode), '-o', output],
                 stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    proc.wait()
    proc = Popen(['bcftools', 'index', '-f', output],
                 stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    proc.wait()
    return proc.poll()


def mean(x):
    if x:
        return round(sum(x) / len(x))


def merge(jasmine_vcfs, records, output):
    all_records = {}
    for recs in records.values():
        for sv_type, rs in recs.items():
            if sv_type == 'HEADER':
                continue
            if sv_type not in all_records:
                all_records[sv_type] = {}
            all_records[sv_type].update(rs)

    merged = []
    for sv_type in jasmine_vcfs:
        vf = VariantFile(jasmine_vcfs[sv_type])
        info_out = '{}/{}.info.bcf'.format(WORK, sv_type)
        vf_out = VcfWriter(info_out, vf.header, remove=REMOVE, drop=True, sort=True)
        for rec in vf:
            IDS = rec.info.get('IDLIST')
            recs = [all_records[sv_type][x] for x in IDS]
            rec.ref = 'N'
            rec.pos = mean([x.pos for x in recs])
            rec.stop = mean([x.stop for x in recs])
            rec.info['SVLEN'] = (mean([x.info['SVLEN'][0] for x in recs if 'SVLEN' in x.info]),)
            rec.filter.clear()
            if any(['PASS' in x.filter for x in recs]):
                rec.filter.add('PASS')
            else:
                for flt in set([y for sl in [x.filter.keys() for x in recs] for y in sl]):
                    rec.filter.add(flt)
            for r in recs:
                r.id = rec.id
                r.pos = rec.pos
                r.ref = 'N'
                r.filter.clear()
                r.filter.add('PASS')
            vf_out.write(rec)
        vf_out.close(index=True)
        sample_outs = []
        for pref in records:
            if sv_type not in records[pref]:
                continue
            pref_out = '{}/{}.{}.clean.bcf'.format(WORK, pref, sv_type)
            header = records[pref]['HEADER']
            vf_out = VcfWriter(pref_out, header, remove=['INFO'], sort=True)
            for rec in records[pref][sv_type].values():
                vf_out.write(rec)
            vf_out.close(index=True)
            sample_outs.append(pref_out)
        out = '{}/{}.merged.bcf'.format(WORK, sv_type)
        cmd = ['bcftools', 'merge', '-m', 'id', '--missing-to-ref', '-Oz', '-o', out, info_out]
        cmd.extend(sample_outs)
        proc = Popen(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.wait()
        proc = Popen(['bcftools', 'index', '-f', out],
                     stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        proc.wait()
        merged.append(out)
    cmd = ['bcftools', 'concat', '-a', '-O{}'.format(MODE[file_ext(output)]), '-o', output]
    cmd.extend(merged)
    proc = Popen(cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    proc.wait()
    return proc.poll()


def main(inputs, output, sv_types, jasmine_dir, jasmine_args):
    jasmine_cp = get_jasmine_cp(jasmine_dir)
    if not os.path.exists(WORK):
        os.makedirs(WORK)
    vcfs, records = proc_vcfs(inputs, sv_types)
    jasmine_vcfs = jasmine(vcfs, jasmine_cp, jasmine_args)
    merge(jasmine_vcfs, records, output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputs', nargs='+', help='input sample vcfs to be merged')
    parser.add_argument('--output', default='output.vcf.gz', help='output variant file')
    parser.add_argument('--work', default='work', help='working directory')
    parser.add_argument('--jasmine-dir', default='/bin', help='path to jasmine jar file')
    parser.add_argument('--sv-type', action='append', help='process only this SVTYPE')
    parser.add_argument('--args', action='append', help='args to pass to jasmine')
    args = parser.parse_args()
    jasmine_args = ['threads=2',
                    'max_dist_linear=0.20',
                    'max_dist=500',
                    'min_overlap=0.80']
    WORK = args.work
    main(args.inputs, args.output, args.sv_type, args.jasmine_dir, jasmine_args)

# TODO: add option --jasmine-args to pass through arbitrary jasmine args
