#!/usr/bin/env python3
"""Staged truvari collapse for a single sample.

Stages:
  1. Collapse high-res callers (Manta, Delly SV, Smoove) among themselves.
  2. Cross-tier collapse: merge stage-1 output with all low-res callers
     (CNVnator, Delly CNV) using loose parameters.
  3. Restore any high-res calls incorrectly dropped in stage 2.

Caller priority is encoded via QUAL scores (100 for highest priority, decreasing
by 1 per caller in the order given), so --keep maxqual preserves the intended
representative.
"""

import argparse
import logging
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s',
                    stream=sys.stderr)
log = logging.getLogger(__name__)


def run(cmd, **kwargs):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    return subprocess.run(cmd, check=True, **kwargs)


def count_records(vcf_path):
    r = run(['bcftools', 'view', '-H', str(vcf_path)], capture_output=True, text=True)
    return len(r.stdout.splitlines())


def set_qual(vcf_path, pass_qual, fail_qual, out_path):
    """Set QUAL by PASS/FAIL status; store original QUAL in ORIGQUAL INFO field."""
    awk_prog = (
        f'BEGIN{{OFS="\\t"}}'
        f'/^#CHROM/{{print "##INFO=<ID=ORIGQUAL,Number=1,Type=Float,'
        f'Description=\\"Original QUAL score\\">"; print; next}}'
        f'/^#/{{print;next}}'
        f'{{origq=$6; $8=($8=="."?"ORIGQUAL="origq:"ORIGQUAL="origq";"$8);'
        f'$6=($7=="PASS")?{pass_qual}:{fail_qual}; print}}'
    )
    p1 = subprocess.Popen(['bcftools', 'view', str(vcf_path)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['awk', awk_prog], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(['bcftools', 'view', '-Oz', '-o', str(out_path)], stdin=p2.stdout)
    p2.stdout.close()
    p3.communicate()
    for p in (p1, p2):
        p.wait()
    if any(p.returncode != 0 for p in (p1, p2, p3)):
        raise RuntimeError(f'set_qual pipeline failed for {vcf_path}')

    run(['bcftools', 'index', '-t', str(out_path)])


def preprocess(name, vcf_path, tmpdir, pass_qual, fail_qual):
    """Sort, rename sample to caller name, assign QUAL by pass/fail status.

    Returns (name, Path) or None if no records remain.
    """
    sorted_  = tmpdir / f'{name}.sorted.vcf.gz'
    altonly  = tmpdir / f'{name}.alt.vcf.gz'
    reheaded = tmpdir / f'{name}.rh.vcf.gz'
    out      = tmpdir / f'{name}.vcf.gz'
    sample_f = tmpdir / f'{name}.sample.txt'

    run(['bcftools', 'sort', '-Oz', '-o', str(sorted_), str(vcf_path)])
    run(['bcftools', 'view', '--min-ac', '1', '-Oz', '-o', str(altonly), str(sorted_)])

    sample_f.write_text(name + '\n')
    run(['bcftools', 'reheader', '-s', str(sample_f), str(altonly), '-o', str(reheaded)])

    set_qual(reheaded, pass_qual, fail_qual, out)

    n = count_records(out)
    if n == 0:
        log.info(f'{name}: no records, skipping')
        return None
    log.info(f'{name}: {n} records (PASS QUAL={pass_qual}, FAIL QUAL={fail_qual})')
    return (name, out)


def merge(vcf_paths, out_path):
    """bcftools merge -m id; if only one VCF, copy it directly."""
    if len(vcf_paths) == 1:
        shutil.copy(str(vcf_paths[0]), str(out_path))
        shutil.copy(str(vcf_paths[0]) + '.tbi', str(out_path) + '.tbi')
    else:
        run(['bcftools', 'merge', '-m', 'id', '-Oz', '-o', str(out_path)]
            + [str(v) for v in vcf_paths])
        run(['bcftools', 'index', '-t', str(out_path)])


def collapse(merged, ref, args_str, out_kept, out_removed):
    """truvari collapse --intra --keep maxqual --chain with caller-supplied args."""
    extra = shlex.split(args_str)
    run([
        'truvari', 'collapse',
        '--intra',
        '--chain',
        '--keep',           'maxqual',
        '--reference',      str(ref),
        '--input',          str(merged),
        '--output',         str(out_kept),
        '--removed-output', str(out_removed),
    ] + extra)
    # truvari writes gzip, not bgzip; recompress so bcftools can index it
    tmp = Path(str(out_kept) + '.tmp')
    run(['bcftools', 'view', '-Oz', '-o', str(tmp), str(out_kept)])
    tmp.replace(out_kept)
    run(['bcftools', 'index', '-t', str(out_kept)])


def get_ids(vcf_path):
    """Return set of non-empty variant IDs from a VCF."""
    result = run(['bcftools', 'query', '-f', '%ID\n', str(vcf_path)],
                 capture_output=True, text=True)
    return {s for s in (line.strip() for line in result.stdout.splitlines())
            if s and s != '.'}


def _popen_awk_pipeline(vcf_path, awk_prog, out_path):
    """bcftools view | awk | bcftools view -Oz; raises on any non-zero return."""
    p1 = subprocess.Popen(['bcftools', 'view', str(vcf_path)], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['awk', awk_prog], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(['bcftools', 'view', '-Oz', '-o', str(out_path)], stdin=p2.stdout)
    p2.stdout.close()
    p3.communicate()
    for p in (p1, p2):
        p.wait()
    if any(p.returncode != 0 for p in (p1, p2, p3)):
        raise RuntimeError(f'awk pipeline failed for {vcf_path}')


def tag_nc1(vcf_path, out_path):
    """Add NC1=<NumCollapsed> to INFO and reset NumCollapsed=0 for stage-2 separation."""
    awk_prog = (
        'BEGIN{OFS="\\t"}'
        '/^#CHROM/{'
        'print "##INFO=<ID=NC1,Number=1,Type=Integer,'
        'Description=\\"Stage-1 NumCollapsed\\">";'
        'print;next}'
        '/^#/{print;next}'
        '{n=split($8,f,";");nc=0;new="";'
        'for(i=1;i<=n;i++){'
        'if(f[i]~/^NumCollapsed=/){v=f[i];sub(/^NumCollapsed=/,"",v);nc=v+0}'
        'else if(f[i]!=".") new=(new==""?f[i]:new";"f[i])}'
        'prefix="NC1="nc";NumCollapsed=0";'
        '$8=(new==""?prefix:prefix";"new);print}'
    )
    _popen_awk_pipeline(vcf_path, awk_prog, out_path)
    run(['bcftools', 'index', '-t', str(out_path)])


def finalize(vcf_path):
    """Restore QUAL; compute N_CALLERS from NC1+NumCollapsed; strip truvari INFO fields."""
    tmp = Path(str(vcf_path) + '.finalize_tmp.vcf.gz')
    awk_prog = (
        'BEGIN{OFS="\\t";drop["NumConsolidated"]=1;drop["CollapseId"]=1}'
        '/^##INFO=<ID=ORIGQUAL/{next}'
        '/^##INFO=<ID=NC1/{next}'
        '/^##INFO=<ID=NumCollapsed/{next}'
        '/^##INFO=<ID=NumConsolidated/{next}'
        '/^##INFO=<ID=CollapseId/{next}'
        '/^#CHROM/{'
        'print "##INFO=<ID=N_CALLERS,Number=1,Type=Integer,'
        'Description=\\"Total caller calls supporting this variant\\">";'
        'print;next}'
        '/^#/{print;next}'
        '{n=split($8,f,";");oq="";nc1=-1;nc2=0;new="";'
        'for(i=1;i<=n;i++){'
        'k=f[i];sub(/=.*/,"",k);'
        'if(f[i]~/^ORIGQUAL=/){oq=f[i];sub(/^ORIGQUAL=/,"",oq)}'
        'else if(f[i]~/^NC1=/){v=f[i];sub(/^NC1=/,"",v);nc1=v+0}'
        'else if(f[i]~/^NumCollapsed=/){v=f[i];sub(/^NumCollapsed=/,"",v);nc2=v+0}'
        'else if(k in drop){}'
        'else if(f[i]!=".") new=(new==""?f[i]:new";"f[i])}'
        'if(oq!="")$6=oq;'
        'total=(nc1>=0)?(nc1+1)+nc2:nc2+1;'
        'nc_str="N_CALLERS="total;'
        '$8=(new==""?nc_str:nc_str";"new);print}'
    )
    _popen_awk_pipeline(vcf_path, awk_prog, tmp)
    tmp.replace(vcf_path)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ref',          required=True)
    p.add_argument('--output',       required=True)
    p.add_argument('--high-res',     nargs='+', required=True, metavar='NAME:VCF')
    p.add_argument('--low-res',      nargs='*', default=[],    metavar='NAME:VCF')
    p.add_argument('--hires-args',
                   default='--refdist 500 --pctseq 0.7 --pctsize 0.7 --pctovl 0.0')
    p.add_argument('--lowres-args',
                   default='--refdist 2000 --pctseq 0.0 --pctsize 0.0 --pctovl 0.8')
    args = p.parse_args()

    def parse_pairs(lst):
        return [item.split(':', 1) for item in lst]

    hires_pairs  = parse_pairs(args.high_res)
    lowres_pairs = parse_pairs(args.low_res or [])

    nh, nl = len(hires_pairs), len(lowres_pairs)
    # QUAL priority: hi-res PASS > hi-res FAIL > lo-res PASS > lo-res FAIL
    # hi-res caller i: pass=100-i, fail=100-nh-i
    # lo-res caller j: pass=100-2*nh-j, fail=100-2*nh-nl-j
    hires_quals  = [(100 - i,      100 - nh - i)      for i in range(nh)]
    lowres_quals = [(100 - 2*nh - j, 100 - 2*nh - nl - j) for j in range(nl)]

    log.info('Caller QUAL scores (PASS / FAIL):')
    for (name, _), (pq, fq) in zip(hires_pairs,  hires_quals):
        log.info(f'  {name} (high-res): {pq} / {fq}')
    for (name, _), (pq, fq) in zip(lowres_pairs, lowres_quals):
        log.info(f'  {name} (low-res):  {pq} / {fq}')

    workdir = Path('.')

    with tempfile.TemporaryDirectory(prefix='preprocess_', dir='.') as _tmp:
        tmpdir = Path(_tmp)

        # --- Preprocess (per-caller files stay in tmpdir, auto-cleaned) ---
        hires_ready  = [r for r in (
            preprocess(n, v, tmpdir, pq, fq)
            for (n, v), (pq, fq) in zip(hires_pairs, hires_quals)
        ) if r]
        lowres_ready = [r for r in (
            preprocess(n, v, tmpdir, pq, fq)
            for (n, v), (pq, fq) in zip(lowres_pairs, lowres_quals)
        ) if r]

        # --- All high-res callers empty: write empty output ---
        if not hires_ready:
            log.info('No high-res records — writing empty output')
            run(['bcftools', 'view', '-h', hires_pairs[0][1], '-Oz', '-o', args.output])
            run(['bcftools', 'index', '-t', args.output])
            return

        # --- Stage 1: high-res intra-collapse (stage files land in workdir) ---
        log.info(f'Stage 1: collapsing {len(hires_ready)} high-res caller(s)')
        s1_merged  = workdir / 'stage1_merged.vcf.gz'
        s1_kept    = workdir / 'stage1_kept.vcf.gz'
        s1_removed = workdir / 'stage1_removed.vcf.gz'
        merge([path for _, path in hires_ready], s1_merged)
        collapse(s1_merged, args.ref, args.hires_args, s1_kept, s1_removed)
        log.info(f'Stage 1 done: {count_records(s1_kept)} kept, '
                 f'{count_records(s1_removed)} collapsed')

        # Tag s1_kept with NC1 (saves NumCollapsed) and reset NumCollapsed=0
        # so stage-2 truvari only counts new cross-tier collapses.
        s1_tagged = workdir / 'stage1_tagged.vcf.gz'
        tag_nc1(s1_kept, s1_tagged)

        # --- No low-res callers: stage 1 is final ---
        if not lowres_ready:
            log.info('No low-res records — stage 1 output is final')
            shutil.copy(str(s1_tagged), args.output)
            finalize(Path(args.output))
            run(['bcftools', 'index', '-t', args.output])
            return

        # --- Stage 2: cross-tier collapse ---
        log.info(f'Stage 2: cross-tier collapse ({len(lowres_ready)} low-res caller(s))')
        s2_merged  = workdir / 'stage2_merged.vcf.gz'
        s2_kept    = workdir / 'stage2_kept.vcf.gz'
        s2_removed = workdir / 'stage2_removed.vcf.gz'
        merge([s1_tagged] + [path for _, path in lowres_ready], s2_merged)
        collapse(s2_merged, args.ref, args.lowres_args, s2_kept, s2_removed)
        log.info(f'Stage 2 done: {count_records(s2_kept)} kept, '
                 f'{count_records(s2_removed)} collapsed')

        # --- Step 3: restore high-res calls incorrectly dropped in stage 2 ---
        dropped = get_ids(s2_removed) & get_ids(s1_kept)
        log.info(f'Step 3: {len(dropped)} high-res ID(s) to restore')

        if dropped:
            ids_file = tmpdir / 'dropped_ids.txt'
            ids_file.write_text('\n'.join(sorted(dropped)))
            dropped_vcf = tmpdir / 'dropped.vcf.gz'
            run(['bcftools', 'view', '-i', f'ID=@{ids_file}',
                 '-Oz', '-o', str(dropped_vcf), str(s1_tagged)])
            run(['bcftools', 'index', '-t', str(dropped_vcf)])
            unsorted = tmpdir / 'unsorted.vcf.gz'
            run(['bcftools', 'concat', '-a', '-Oz', '-o', str(unsorted),
                 str(s2_kept), str(dropped_vcf)])
            run(['bcftools', 'sort', '-T', '.', '-Oz', '-o', args.output, str(unsorted)])
        else:
            shutil.copy(str(s2_kept), args.output)

    finalize(Path(args.output))
    run(['bcftools', 'index', '-t', args.output])
    log.info(f'Output: {count_records(args.output)} variants → {args.output}')


if __name__ == '__main__':
    main()
