#!/usr/bin/env python3
"""End-to-end synthetic regression tests for v11.

Requires numpy/pandas/Biopython; substitutes only pyfaidx in subprocess tests
because pyfaidx may be unavailable in a development sandbox.
No production script uses the stub.
"""
from __future__ import annotations
import csv
import importlib.util
import sys
import types
from pathlib import Path
import tempfile
import test_04_build_compleasm_feature_tsvs_v10 as base

SCRIPT=Path(__file__).resolve().with_name('04_build_compleasm_feature_tsvs_v11.py')
NULL=base.NULL
ODB='414176at8457'
ACC='GWHFILM00000000.1'
CDS='ATGCCGTTAAAC'


def fixture(tmp, *, chosen=0, first_invalid=False, same_scaffold=False,
            same_sequence=False, positive=False, with_log=True,
            conflict=False, tsv_mismatch=False, genome_mismatch=False,
            multiple_tsv=False, unmapped=False, missing_qc=False):
    v=base.prepare(tmp, qc_dup=True)
    origacc=v['acc']
    # Replace IDs in authoritative identity files, not the old unrelated
    # genomic sequences. The test genome is a compact reference assembly.
    for path in (v['manifest'], v['genome_table'], v['qc_path'],
            v['root']/'records/compleasm/records/metadata.csv',
            v['root']/'records/genomes_metadata.csv'):
        path.write_text(path.read_text().replace(origacc, ACC).replace('32015at8457', ODB))
    (Path(tmp)/'orthodb'/'links_to_ODB12.txt').write_text(
        f'111at8457\tOrtholog 111\turl\n{ODB}\tOrtholog {ODB}\turl\n')
    source=v['a']
    full=v['f1']
    exons = '5_10_-|21_26_-' if not positive else '5_10_+|21_26_+'
    strand = '+' if positive else '-'
    if positive:
        chrD='T'*4+'ATGCCG'+'A'*10+'TTAAAC'+'A'*40
    else:
        chrD='T'*4+'GTTTAA'+'A'*10+'CGGCAT'+'A'*40
    if genome_mismatch:
        chrD=chrD[:4]+'X'+chrD[5:]
    chrE='G'*90
    chrF='C'*90
    # Upstream copy can be malformed, so the extractor moves to a later one.
    rows=[]
    for i,chrom in enumerate(['chrD','chrE','chrF']):
        if first_invalid and i==0:
            row=[ODB,'Duplicated','missing_chr',23000,12,strand,5,26,exons,0.9]
        else:
            row=[ODB,'Duplicated',chrom,23000+i,12,strand,5,26,exons,0.9]
        rows.append(row)
    if first_invalid:
        rows[1]=[ODB,'Duplicated','chrD',23001,12,strand,5,26,exons,0.9]
        rows[2]=[ODB,'Duplicated','chrF',23002,12,strand,5,26,exons,0.9]
    if same_scaffold:
        # Create a second *distinct* locus on chrD that produces same CDS.
        chrD += 'A'*18 + (('ATGCCG'+'A'*10+'TTAAAC') if positive else ('GTTTAA'+'A'*10+'CGGCAT'))
        # Original length=66, +18 -> second exons start 85 and 101.
        rows[1]=[ODB,'Duplicated','chrD',23001,12,strand,85,106,
                 ('85_90_+|101_106_+' if positive else '85_90_-|101_106_-'),0.9]
        chosen=0
    if same_sequence and not same_scaffold:
        chrE=chrD
    header=['Gene','Status','Sequence','Score','Length','Strand','Gene Start','Gene End','Codons','Identity']
    allrows=[['111at8457','Single','chr1',500,12,'+',10,45,'10_21_+|30_45_+',0.9],*rows]
    base.write_table(full,header,allrows)
    v['g1'].write_text('>chr1\n'+'GCT'*200+'\n>chrD\n'+chrD+'\n>chrE\n'+chrE+'\n>chrF\n'+chrF+'\n')
    record=rows[chosen]
    loc=record[2]
    new_cds=source/'Indotyphlops_braminus_cds_compleasm.fasta'
    new_cds.write_text(f'>111at8457\nGCTGCTGCTGCT\n>{ODB}\n{CDS}\n')
    tsv=new_cds.with_suffix('.tsv')
    tsvrows=[['111at8457','indotypholops_braminus','chr1','+','GCTGCTGCTGCT'],
             [ODB,'indotypholops_braminus',loc,strand,CDS if not tsv_mismatch else 'T'+CDS[1:]]]
    if multiple_tsv: tsvrows.append(tsvrows[-1])
    base.write_table(tsv,['Busco ID','species','location','strand','sequence'],tsvrows)
    log=new_cds.with_name(new_cds.name.replace('.fasta','_records.txt'))
    if with_log:
        log.write_text(f'Wrote {ODB} (Duplicated) from '
                       f'{"wrong_chrom" if conflict else loc}, strand {strand}, length {len(CDS)} bp.\n'
                       f'Skipping additional Duplicated copy of {ODB}: a valid sequence for this BUSCO ID has already been written.\n')
    meta=v['root']/'records/compleasm/records/metadata.csv'
    meta.write_text(meta.read_text().replace(str(v['c1']),str(new_cds)))
    v['sequence_table'].write_text(v['sequence_table'].read_text().replace('QEPC01013744.1','chrD'))
    if not unmapped:
        with v['sequence_table'].open('a') as h:
            h.write('4\t1\tchrE\n5\t1\tchrF\n')
    else:
        v['sequence_table'].write_text(v['sequence_table'].read_text().replace('chrD','not_chrD'))
    if missing_qc:
        qrows=base.read_rows(v['qc_path'])
        qcols=list(qrows[0])
        base.write_table(v['qc_path'],qcols,[[r[c] for c in qcols] for r in qrows if r['odb12_id']!=ODB])
    v.update(selected=rows[chosen], cds_tsv=tsv, cds_fasta=new_cds,
             extraction_log=log, full=full)
    return v


def execute(info, run=False):
    return base.run(info,run=run,script=SCRIPT)


def rows(info,name):
    return base.read_rows(info['out']/name)


def successful(info,why):
    p=execute(info)
    base.assert_success(p,why)
    orth=[r for r in rows(info,'orthologs.tsv') if r['odb12_id']==ODB]
    assert len(orth)==1
    qc=[r for r in rows(info,'ortholog_initial_qc.tsv') if r['ortholog_pk']==orth[0]['ortholog_pk']]
    assert len(qc)==1 and qc[0]['genome_pk']==orth[0]['genome_pk']
    return p,orth[0]


def test_observed_source_geometry_synthetic():
    """Exercise documented lines/coords/13950-bp length with virtual genome.

    The *real* 85-genome FASTA/full_table were not uploaded. This uses the
    documented coordinates and source identifiers, with synthetic exon bases.
    """
    if 'pyfaidx' not in sys.modules:
        fake=types.ModuleType('pyfaidx')
        fake.Fasta=object
        fake.FaidxException=type('FaidxException',(Exception,),{})
        sys.modules['pyfaidx']=fake
    spec=importlib.util.spec_from_file_location('compleasm_v11_source_test',SCRIPT)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=mod
    spec.loader.exec_module(mod)
    class VirtualChrom:
        def __init__(self,base): self.base=base
        def __len__(self): return 200_000_000
        def __getitem__(self,sl): return self.base*(sl.stop-sl.start)
    class VirtualGenome:
        chroms={'GWHFILM00000003.1': VirtualChrom('A'),
                'GWHFILM00000001.1': VirtualChrom('C'),
                'GWHFILM00000002.1': VirtualChrom('G')}
        def __getitem__(self,x): return self.chroms[x]
    with tempfile.TemporaryDirectory() as d:
        root=Path(d)
        odb='414176at8457'
        seq='T'*13950
        fasta=root/'Indotyphlops_braminus_cds_compleasm.fasta'
        fasta.write_text(f'>{odb}\n{seq}\n')
        base.write_table(fasta.with_suffix('.tsv'),
           ['Busco ID','species','location','strand','sequence'],
           [[odb,'indotypholops_braminus','GWHFILM00000003.1','-',seq]])
        fasta.with_name(fasta.stem+'_records.txt').write_text(
           f'Wrote {odb} (Duplicated) from GWHFILM00000003.1, strand -, length 13950 bp.\n')
        # This is a synthetic two-exon geometry spanning the documented gene
        # boundaries, not a claim about the real exon structure.
        candidate_rows=[]
        for ln,chrom,lo,hi,score in [
            (3,'GWHFILM00000003.1',174490142,174670948,23130),
            (4,'GWHFILM00000001.1',177078927,177263988,23157),
            (5,'GWHFILM00000002.1',176145762,176329170,23107)]:
            codons=f'{lo}_{lo+6974}_-|{hi-6974}_{hi}_-'
            record=mod.FullRecord(odb,'Duplicated',chrom,'-',lo,hi,[])
            candidate_rows.append((record,{'source_line':ln,
                'source_sequence_id':chrom,'source_strand':'-',
                'source_codon_intervals':codons}))
        ctx=mod.load_selected_copy_provenance(fasta)
        idx,evidence=mod.resolve_historical_copy(ACC,odb,candidate_rows,ctx,VirtualGenome(),{odb:seq})
        assert idx==0 and evidence['selected_source_line']==3
        assert evidence['reconstructed_cds_length']==13950
        r=candidate_rows[idx][0]
        assert (r.sequence_id,r.start,r.end,r.strand)==(
            'GWHFILM00000003.1',174490142,174670948,'-')
        print('PASS synthetic observed coordinates: source 3, chr3, 174490142..174670948, negative strand, 13950 bp')


def test_verified_default():
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d)
        p, r=successful(x,'negative strand triploid')
        assert (r['sequence_pk'],r['strand'],r['start'],r['end'],r['status'])==('2','-','5','26','Duplicated'),r
        assert r['gc']!=NULL and r['gc3']!=NULL and r['gc4']!=NULL
        da=[row for row in rows(x,'ortholog_duplicate_full_table_audit.tsv') if row['odb12_id']==ODB]
        assert len(da)==3
        assert [r['genomic_copy_selected'] for r in da]==['TRUE','FALSE','FALSE']
        assert [r['source_line'] for r in da]==['3','4','5']
        assert all(r['resolution_status']=='eligible_duplicate_provenance_verified' for r in da)
        assert da[0]['selected_source_line']=='3' and da[0]['cds_sequence_match']=='TRUE'
        ai=[r for r in rows(x,'ortholog_source_audit.tsv') if r['odb12_id']==ODB][0]
        assert ai['mapping_status']=='mapped_genomic_sequence' and 'verified' in ai['mapping_notes']
        assert any(intron['ortholog_pk']==r['ortholog_pk'] for intron in rows(x,'intron_compleasm.tsv'))
        assert 'Provenance-verified Indotyphlops groups: 1' in p.stdout
        print('PASS: provenance-verified negative strand, source line 3, 1:1 SQL, intron and flank')


def test_later_copy():
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d, chosen=1, first_invalid=True)
        _,r=successful(x,'fallback after unsuccessful candidate')
        assert r['sequence_pk']=='2'
        da=[row for row in rows(x,'ortholog_duplicate_full_table_audit.tsv') if row['odb12_id']==ODB]
        assert [r['genomic_copy_selected'] for r in da]==['FALSE','TRUE','FALSE']
        assert da[1]['selected_source_line']=='4'
        print('PASS: earlier extraction fails; later verified copy has original source line 4')


def test_positive():
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d, positive=True)
        _,r=successful(x,'positive strand')
        assert r['strand']=='+'
        print('PASS: positive-strand reconstruction')


def test_no_log():
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d,with_log=False)
        _,r=successful(x,'TSV + FASTA without optional log')
        assert r['sequence_pk']=='2'
        print('PASS: absent optional log does not undermine three-source sequence verification')


def test_unmapped():
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d,unmapped=True)
        _,r=successful(x,'unmapped verified source')
        assert all(r[c]==NULL for c in ('sequence_pk','start','end','strand'))
        assert all(r[c]!=NULL for c in ('gc','gc3','gc4'))
        src=[r for r in rows(x,'ortholog_source_audit.tsv') if r['odb12_id']==ODB][0]
        assert 'verified' in src['mapping_notes']
        print('PASS: verified but unmapped selected source preserves CDS composition')


def expect_failure(**kw):
    with tempfile.TemporaryDirectory() as d:
        x=fixture(d,**kw)
        x['out'].mkdir()
        sentinel=x['out']/'sentinel.txt'
        sentinel.write_text('prior output unchanged')
        p=execute(x)
        assert p.returncode!=0 and p.stderr,(kw,p.stdout[-1000:],p.stderr[-1000:])
        assert sentinel.read_text()=='prior output unchanged'
        assert not (x['out']/'orthologs.tsv').exists()
        assert (Path(d)/'output.v11_failure_report.tsv').exists(),p.stderr
        failure=base.read_rows(Path(d)/'output.v11_duplicate_failure_audit.tsv')
        assert len(failure)>=3
        return p.stderr, [r['resolution_status'] for r in failure if r['odb12_id']==ODB]


def test_failures():
    for label,params,code in [
        ('same-scaffold copied CDS',dict(same_scaffold=True),'eligible_duplicate_ambiguous_source'),
        ('contradictory log',dict(conflict=True),'eligible_duplicate_conflicting_provenance'),
        ('mismatched extraction TSV sequence',dict(tsv_mismatch=True),'eligible_duplicate_cds_mismatch'),
        ('mismatched original genomic sequence',dict(genome_mismatch=True),'eligible_duplicate_cds_mismatch'),
        ('historical TSV ambiguous duplicate identifier',dict(multiple_tsv=True),'eligible_duplicate_conflicting_provenance'),
        ('historical first-successful-copy conflict',dict(chosen=1,same_sequence=True),'eligible_duplicate_conflicting_provenance'),
        ('missing required initial QC',dict(missing_qc=True),'QC-eligible'),
    ]:
        stderr,codes=expect_failure(**params)
        assert code in stderr,(label,stderr[-1700:])
        if code.startswith('eligible_duplicate_'):
            assert codes and all(x==code for x in codes), (label,codes)
        print('PASS failure:',label)


def test_unambiguous_v10_same_outputs():
    with tempfile.TemporaryDirectory() as d:
        x=base.prepare(d, include_second=False)
        out10=Path(d)/'out_v10'
        out11=Path(d)/'out_v11'
        x['args'][-1]=str(out10)
        p10=base.run(x,script=SCRIPT.with_name('04_build_compleasm_feature_tsvs_v10.py'))
        base.assert_success(p10,'v10 reference')
        outputs10={p.name:base.read_rows(p) for p in out10.glob('*.tsv')}
        x['args'][-1]=str(out11)
        p11=execute(x)
        base.assert_success(p11,'v11 reference')
        outputs11={p.name:base.read_rows(p) for p in out11.glob('*.tsv')}
        assert outputs10.keys()==outputs11.keys()
        for filename in outputs10:
            if filename=='ortholog_duplicate_full_table_audit.tsv': continue
            a,b=outputs10[filename],outputs11[filename]
            for row in a+b:
                row.pop('created_at',None)
            assert a==b,filename
        print('PASS: all v10 SQL and feature data unchanged for unambiguous sources')


def test_unresolved_and_run_qc():
    with tempfile.TemporaryDirectory() as d:
        x=base.prepare(d)
        base.setup_run(x)
        p=execute(x,run=True)
        base.assert_success(p,'v10 run QC preservation')
        audit=rows(x,'ortholog_duplicate_full_table_audit.tsv')
        assert len(audit)==2 and all(r['resolution_status']=='ambiguous_ineligible_retained' for r in audit)
        run=rows(x,'ortholog_run_qc.tsv')
        assert len(run)==2 and all(r['included_in_nhphyml_alignment']=='TRUE' for r in run)
        print('PASS: QC-ineligible duplicates and existing configured run QC retained')


if __name__=='__main__':
    test_observed_source_geometry_synthetic()
    test_verified_default()
    test_later_copy()
    test_positive()
    test_no_log()
    test_unmapped()
    test_failures()
    test_unambiguous_v10_same_outputs()
    test_unresolved_and_run_qc()
    print('ALL V11 SYNTHETIC TESTS PASSED')
