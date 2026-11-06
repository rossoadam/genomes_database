"""Local regression harness; uses test-only pyfaidx stub if package unavailable.
Never distributed as production runtime dependency."""
from __future__ import annotations
import csv
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT_SCRIPT=Path(__file__).resolve().with_name('04_build_compleasm_feature_tsvs_v10.py')
V9_SCRIPT=ROOT_SCRIPT.with_name('04_build_compleasm_feature_tsvs_v9.py')
NULL=r'\N'


def write_table(path, cols, rows, delimiter='\t'):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as h:
        w=csv.writer(h,delimiter=delimiter)
        w.writerow(cols)
        w.writerows(rows)


def read_rows(path):
    with path.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))


def prepare(base, *, sphenodon_status='Duplicated', second_status='Duplicated', identical=False, qc_dup=False, single_unmapped=False, include_second=True):
    root=Path(base)/'genomes'
    root.mkdir()
    acc='GCA_003113815.1'
    other='GCA_999999999.1'
    a= root/'records'/ 'compleasm' /'sources'; a.mkdir(parents=True)
    # pyfaidx dependency is missing from the test environment; inject a
    # faithful tiny FASTA subset implementation only into the test process.
    stub=Path(base)/'stubs';stub.mkdir()
    (stub/'pyfaidx.py').write_text('''class FaidxException(Exception): pass
class Fasta:
    def __init__(self, path, as_raw=True, build_index=True):
        self.seq={}
        key=None
        for line in open(path):
            line=line.strip()
            if line.startswith('>'):
                key=line[1:].split()[0]
                self.seq[key]=''
            elif key: self.seq[key]+=line
    def __contains__(self, key):return key in self.seq
    def __getitem__(self, key):return self.seq[key]
''')
    full_cols=['Gene','Status','Sequence','Score','Length','Strand','Gene Start','Gene End','Codons','Identity']
    full=[['111at8457','Single', 'chr1','510','898','+',10,45,'10_21_+|30_45_+','99.5'],
          ['32015at8457',sphenodon_status,'QEPC01013744.1','4260','898','+',268123,362175,'268123_268500_+','98'],
          ['32015at8457',second_status,'QEPC01013744.1','4075','898','+',387400,510339,'387400_387800_+','97']]
    if identical: full[-1]=full[-2][:]
    if not include_second: full=full[:-1]
    f1=a/'sphenodon_full.tsv';write_table(f1,full_cols,full)
    f2=a/'other_full.tsv';write_table(f2,full_cols,[['111at8457','Single','chr1', '500','898','+',30,60,'30_36_+|44_60_+','99']])
    c1=a/'sphenodon_cds.fasta';c1.write_text('>111at8457\nGCTGCTGCTGCT\n>32015at8457\nGCTGCCGCAGCT\n')
    c2=a/'other_cds.fasta';c2.write_text('>111at8457\nGCTGCCGCTGCT\n')
    g1=a/'sphenodon_genome.fna';g1.write_text('>chr1\n'+'GCT'*200+'\n>QEPC01013744.1\n'+'GCT'*1000+'\n')
    g2=a/'other_genome.fna';g2.write_text('>chr1\n'+'GCT'*200+'\n')
    write_table(root/'records'/'compleasm'/'records'/'metadata.csv', ['accession','cds_fasta','full_table'],
                [[acc,str(c1),str(f1)],[other,str(c2),str(f2)]],delimiter=',')
    write_table(root/'records'/'genomes_metadata.csv',['accession','path_to_fna'],
                [[acc,str(g1)],[other,str(g2)]],delimiter=',')
    manifest=root/'records'/'project_manifests'/'manifest.csv';write_table(manifest,['accession','species_key'],[[acc,'Sphenodon_punctatus'],[other,'Test_species']], delimiter=',')
    genome_table=Path(base)/'genomes.tsv';write_table(genome_table,['genome_pk','accession_id','species_pk'],[[1,acc,1],[2,other,2]])
    sequence_table=Path(base)/'sequences.tsv';write_table(sequence_table,['sequence_pk','genome_pk','sequence_id'],[[1,1,'chr1'],[2,1,'QEPC01013744.1'],[3,2,'chr1' if not single_unmapped else 'other_chr']])
    qc_cols=['accession','odb12_id','final_qc','macse_qc','raw_cds_qc','raw_length','raw_length_mod_3','raw_internal_stop_count','raw_terminal_stop','raw_ambiguous_fraction','raw_invalid_fraction','raw_failure_reasons','macse_frameshift','macse_internal_stop','macse_terminal_stop','macse_failure_reasons','qc_version']
    def qcrow(accession, odb):
        return [accession,odb,NULL,NULL,'pass',12,0,0,'FALSE',0,0,NULL,NULL,NULL,NULL,NULL,'v4']
    qrows=[qcrow(acc,'111at8457'),qcrow(other,'111at8457')]
    if qc_dup:qrows.append(qcrow(acc,'32015at8457'))
    qc_path=Path(base)/'ortholog_validity.tsv';write_table(qc_path,qc_cols,qrows)
    orthodb=Path(base)/'orthodb';orthodb.mkdir();(orthodb/'links_to_ODB12.txt').write_text('111at8457\tOrtholog 111\turl\n32015at8457\tOrtholog 32015\turl\n')
    out=Path(base)/'output';
    args=['--genomes',str(root),'--manifest',str(manifest),'--orthodb',str(orthodb),
          '--sequences-tsv',str(sequence_table),'--genomes-tsv',str(genome_table),
          '--ortholog-validation',str(qc_path),'--outdir',str(out)]
    return locals()


def run(info, *, run=False, script=ROOT_SCRIPT):
    env=os.environ.copy();env['PYTHONPATH']=str(info['stub'])+os.pathsep+env.get('PYTHONPATH','')
    args=[sys.executable,str(script),*info['args']]
    if run:
        args+=['--run-config',str(info['run_cfg']),'--nhphyml-runs',str(info['runs_tsv'])]
    return subprocess.run(args, env=env,capture_output=True,text=True,timeout=50)


def setup_run(info, *, incl_dup=False, eval_dup=False):
    base=Path(info['base']);inp=base/'nhphyml_inputs';inp.mkdir(exist_ok=True)
    inp.joinpath('111at8457.fasta').write_text('>Sphenodon_punctatus\nGCTGCTGCTGCT\n>Test_species\nGCTGCCGCTGCT\n')
    if not incl_dup and inp.joinpath('32015at8457.fasta').exists():
        inp.joinpath('32015at8457.fasta').unlink()
    if incl_dup:
        inp.joinpath('32015at8457.fasta').write_text('>Sphenodon_punctatus\nGCTGCTGCTGCT\n')
    ia=base/'included_accessions.txt';ia.write_text(info['acc']+'\n'+info['other']+'\n')
    so=base/'shared_orthologs.txt';so.write_text('111at8457\n'+('32015at8457\n' if incl_dup else ''))
    vali=base/'final_macse_validation.tsv';rows=[[info['acc'],'111at8457','pass',NULL,'vfinal'],[info['other'],'111at8457','pass',NULL,'vfinal']]
    if eval_dup:rows.append([info['acc'],'32015at8457','fail','frameshift','vfinal'])
    write_table(vali,['accession','odb12_id','final_alignment_qc','macse_failure_reasons','qc_version'],rows)
    info['runs_tsv']=base/'nhphyml_runs.tsv';write_table(info['runs_tsv'],['nhphyml_run_pk','run_name','alignment_name'],[[1,'80','80_2026oct14']])
    info['run_cfg']=base/'run_config.tsv';write_table(info['run_cfg'],['nhphyml_run_pk','run_name','final_macse_validation','included_accessions','shared_orthologs','nhphyml_input_dir','qc_version'],[[1,'80',str(vali),str(ia),str(so),str(inp),'vfinal']])


def assert_success(proc, name):
    if proc.returncode != 0:
        raise AssertionError(f'{name} failed:\n{proc.stdout[-2800:]}\n{proc.stderr[-2800:]}')


def assert_failure(proc, pattern, name):
    if proc.returncode==0 or pattern not in proc.stderr:
        raise AssertionError(f'{name} expected failure {pattern}:\n{proc.stdout[-2500:]}\n{proc.stderr[-2500:]}')


def test_sphenodon():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d)
        p=run(info); assert_success(p,'Sphenodon duplicate')
        orth=read_rows(info['out']/'orthologs.tsv')
        dup=[r for r in orth if r['odb12_id']=='32015at8457']
        assert len(dup)==1 and dup[0]['genome_pk']=='1' and dup[0]['status']=='Duplicated'
        for field in ['sequence_pk','strand','start','end','gc','gc3','gc4']:
            assert dup[0][field]==NULL, (field,dup)
        qc=[r for r in read_rows(info['out']/'ortholog_initial_qc.tsv') if r['ortholog_pk']==dup[0]['ortholog_pk']]
        assert len(qc)==1 and all(qc[0][c]==NULL for c in info['qc_cols'][2:])
        audit=[r for r in read_rows(info['out']/'ortholog_duplicate_full_table_audit.tsv') if r['odb12_id']=='32015at8457']
        assert len(audit)==2 and all(r['resolution_status']=='ambiguous_ineligible_retained' and r['genomic_copy_selected']=='FALSE' for r in audit)
        assert [r['source_line'] for r in audit]==['3','4']
        source=[r for r in read_rows(info['out']/'ortholog_source_audit.tsv') if r['odb12_id']=='32015at8457']
        assert len(source)==1 and source[0]['mapping_status']=='ambiguous_duplicated_ineligible'
        assert all(source[0][c]==NULL for c in ['source_sequence_id','source_strand','source_start','source_end','source_codon_intervals'])
        assert not any(r['ortholog_pk']==dup[0]['ortholog_pk'] for r in read_rows(info['out']/'intron_compleasm.tsv'))
        assert not any(r['ortholog_pk']==dup[0]['ortholog_pk'] for r in read_rows(info['out']/'flanks_compleasm.tsv'))
        assert [r for r in read_rows(info['out']/'orthologs.tsv') if r['odb12_id']=='111at8457'][0]['gc'] != NULL
        assert 'SQL NULL serialization validation passed' in p.stdout
        print('PASS A: Sphenodon unresolved duplicate, 1:1 QC, 2-line audit, no derived features')
        return info


def test_identical():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d,identical=True)
        p=run(info);assert_success(p,'identical')
        assert len([r for r in read_rows(info['out']/'orthologs.tsv') if r['odb12_id']=='32015at8457'])==1
        assert all(r['resolution_status']=='identical_records_collapsed' for r in read_rows(info['out']/'ortholog_duplicate_full_table_audit.tsv'))
        print('PASS C: identical repeats collapsed')


def test_fatal_groups():
    for name,kw,expected in [
        ('distinct Single',dict(sphenodon_status='Single',second_status='Single'),'ambiguous_eligible_error'),
        ('mixed status',dict(sphenodon_status='Single',second_status='Duplicated'),'conflicting_status_error'),
        ('contradictory QC',dict(qc_dup=True),'contradictory_qc_error')]:
        with tempfile.TemporaryDirectory() as d:
            info=prepare(d,**kw)
            info['out'].mkdir();(info['out']/'sentinel.txt').write_text('original')
            proc=run(info);assert_failure(proc,expected,name)
            assert (info['out']/'sentinel.txt').read_text()=='original'
            assert not (info['out']/'orthologs.tsv').exists()
            r=read_rows(Path(d)/'output.v10_duplicate_failure_audit.tsv')
            assert len(r)==2 and all(x['resolution_status']==expected for x in r)
            print('PASS',name,': preserved prior output and detailed failure-copy report')


def test_briminus():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d)
        # make both source rows from the eligible triploid accession
        acc='GWHFILM00000000.1'
        for path in [info['manifest'], info['genome_table'],info['qc_path'],info['root']/'records/compleasm/records/metadata.csv',info['root']/'records/genomes_metadata.csv']:
            v=path.read_text().replace(info['acc'],acc);path.write_text(v)
        proc=run(info)
        assert_failure(proc,'ambiguous_eligible_error','Indotyphlops')
        print('PASS E: triploid exception does not select one of two distinct copies')


def test_unmapped_single():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d,single_unmapped=True)
        proc=run(info);assert_success(proc,'unmapped Single')
        ortho=[r for r in read_rows(info['out']/'orthologs.tsv') if r['genome_pk']=='2' and r['odb12_id']=='111at8457'][0]
        assert ortho['sequence_pk']==NULL and ortho['gc']!=NULL and ortho['gc3']!=NULL
        source=[r for r in read_rows(info['out']/'ortholog_source_audit.tsv') if r['genome_pk']=='2' and r['odb12_id']=='111at8457'][0]
        assert source['mapping_status']=='unmapped_sequence_id'
        print('PASS H: unmapped eligible Single retains GC/GC3 and existing mapping classification')


def test_run():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d)
        setup_run(info)
        proc=run(info,run=True);assert_success(proc,'run QC')
        q=read_rows(info['out']/'ortholog_run_qc.tsv')
        assert len(q)==2 and [r['ortholog_run_qc_pk'] for r in q]==['1','2']
        assert [r['included_in_nhphyml_alignment'] for r in q]==['TRUE','TRUE']
        assert len(read_rows(info['out']/'ortholog_run_selection_audit.tsv'))==3
        print('PASS I: per-sequence run QC and audit include excluded duplicate without failing')
        # run validator forbids explicitly evaluated duplicate even if its QC says fail
        setup_run(info,eval_dup=True)
        proc=run(info,run=True)
        assert_failure(proc,'unresolved genomic copies contradict run evidence','run QC contradiction')
        print('PASS I.1: final-MACSE record for excluded duplicate rejected')
        # input inclusion alone rejected
        setup_run(info,incl_dup=True,eval_dup=False)
        proc=run(info,run=True)
        assert_failure(proc,'unresolved genomic copies contradict run evidence','included duplicate')
        print('PASS I.2: nhPhyML alignment inclusion for excluded duplicate rejected')


def test_unchanged_v9_behavior():
    if not V9_SCRIPT.exists():
        print('SKIP v9 side-by-side comparison: place original v9 script beside this test file')
        return
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d,include_second=False)
        info['args'][-1]=str(Path(d)/'out_v9')
        proc9=run(info, script=V9_SCRIPT); assert_success(proc9,'unambiguous v9 baseline')
        info['args'][-1]=str(Path(d)/'out_v10')
        proc10=run(info); assert_success(proc10,'unambiguous v10')
        p9,p10=Path(d)/'out_v9',Path(d)/'out_v10'
        files9={x.name for x in p9.glob('*.tsv')}
        files10={x.name for x in p10.glob('*.tsv')}
        assert files10-files9=={'ortholog_duplicate_full_table_audit.tsv'}
        for name in files9:
            a=read_rows(p9/name)
            b=read_rows(p10/name)
            for row in a+b:
                if 'created_at' in row: row['created_at']='TIMESTAMP'
            assert a==b, f'{name}: v9-v10 differences unrelated to repeated records: {a[:2]} vs {b[:2]}'
        assert not read_rows(p10/'ortholog_duplicate_full_table_audit.tsv')
        print('PASS B: v9/v10 table-by-table equality for unchanged single-copy inputs (excluding creation timestamps)')


def test_other_ineligible_and_bad_status():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d,sphenodon_status='Fragmented',second_status='Fragmented')
        proc=run(info);assert_success(proc,'Fragmented multi-copy')
        row=[x for x in read_rows(info['out']/'orthologs.tsv') if x['odb12_id']=='32015at8457'][0]
        assert row['status']=='Fragmented' and row['sequence_pk']==NULL and row['gc']==NULL
        mapped=[x for x in read_rows(info['out']/'ortholog_source_audit.tsv') if x['odb12_id']=='32015at8457'][0]
        assert mapped['mapping_status']=='ambiguous_ineligible_genomic_copies'
        print('PASS: distinct same-status QC-ineligible Fragmented records unresolved')
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d,sphenodon_status='Unknown',second_status='Unknown')
        proc=run(info);assert_failure(proc,'invalid_status_error','invalid Compleasm status')
        print('PASS: invalid repeated Compleasm status fails closed')


def test_v13_headers():
    orth=['ortholog_pk','genome_pk','sequence_pk','odb12_id','status','strand','start','end','gc','gc3','gc4']
    iq=['ortholog_pk','genome_pk','sequence_pk','odb12_id','final_qc','macse_qc','raw_cds_qc','raw_length','raw_length_mod_3','raw_internal_stop_count','raw_terminal_stop','raw_ambiguous_fraction','raw_invalid_fraction','raw_failure_reasons','macse_frameshift','macse_internal_stop','macse_terminal_stop','macse_failure_reasons','qc_version']
    rq=['ortholog_run_qc_pk','ortholog_pk','nhphyml_run_pk','final_alignment_qc','macse_failure_reasons','included_in_nhphyml_alignment','qc_version']
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d)
        proc=run(info);assert_success(proc,'v13 headers')
        for name,expected in [('orthologs.tsv',orth),('ortholog_initial_qc.tsv',iq),('ortholog_run_qc.tsv',rq)]:
            with (info['out']/name).open() as f:
                observed=next(csv.reader(f,delimiter='\t'))
            assert observed==expected,(name,observed,expected)
        assert len(read_rows(info['out']/'orthologs.tsv'))==len(read_rows(info['out']/'ortholog_initial_qc.tsv'))
        print('PASS J: three v13 SQL schemas, exact header order, one-to-one initial QC')


def test_atomic_replace_on_success():
    with tempfile.TemporaryDirectory() as d:
        info=prepare(d)
        before=run(info);assert_success(before,'initial build')
        first=(info['out']/'orthologs.tsv').read_bytes()
        second=run(info);assert_success(second,'repeated build')
        assert first==(info['out']/'orthologs.tsv').read_bytes()
        assert not list(info['out'].parent.glob('*.v10_previous_build'))
        assert not list(info['out'].parent.glob('.compleasm_features_v10_stage_*'))
        print('PASS: repeated successful build preserves deterministic ortholog rows, cleans up staged/backup folders')


if __name__=='__main__':
    test_sphenodon()
    test_identical()
    test_fatal_groups()
    test_briminus()
    test_unmapped_single()
    test_run()
    test_unchanged_v9_behavior()
    test_other_ineligible_and_bad_status()
    test_v13_headers()
    test_atomic_replace_on_success()
    print('ALL LOCAL V10 TESTS PASSED')
