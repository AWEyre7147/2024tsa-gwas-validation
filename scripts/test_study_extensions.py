"""Meaningful edge-case checks for diagnostic, inventory and power extensions."""
import unittest
import numpy as np
from additional_qc import eligible
from build_gene_inventory import gene_class,attrs,resolve_transcript_genes
from power_framework import benchmark_power,wilson


class ExtensionChecks(unittest.TestCase):
    def test_precursor_and_mature_RNA_gene_relationships(self):
        import sqlite3,json
        db=sqlite3.connect(':memory:')
        db.executescript('CREATE TABLE genes(record_key TEXT); CREATE TABLE transcripts(record_key TEXT,parent_gene_keys_json TEXT);')
        db.execute('INSERT INTO genes VALUES (?)',('gene:1',))
        db.executemany('INSERT INTO transcripts VALUES (?,?)',[('precursor',json.dumps(['gene:1'])),('mature',json.dumps(['precursor']))])
        self.assertEqual(resolve_transcript_genes(db),1)
        self.assertEqual(json.loads(db.execute('SELECT gene_keys_json FROM transcripts WHERE record_key="mature"').fetchone()[0]),['gene:1'])
        db.close()

    def test_structural_selection_excludes_duplicate_coordinates(self):
        data=[['1','a','0','10','A','G'],['1','b','0','10','T','C'],
              ['1','c','0','20','A','C'],['39','x','0','9000000','A','G'],['2','d','0','0','A','T']]
        self.assertEqual(eligible(data),[2])

    def test_distinct_noncoding_pseudogene_and_immune_classes(self):
        self.assertEqual(gene_class('protein_coding'),'protein_coding')
        self.assertEqual(gene_class('lncRNA'),'noncoding')
        self.assertEqual(gene_class('processed_pseudogene'),'pseudogene')
        self.assertEqual(gene_class('IG_V_gene'),'immune_gene_segment')
        self.assertEqual(gene_class('unknown'),'other_or_unclassified')
        self.assertEqual(attrs('ID=gene-x;description=some%3Bdescription','RefSeq')['description'],'some;description')

    def test_null_benchmark_and_effect_direction_symmetry(self):
        p,_=benchmark_power(.2,1,7,39,.05)
        self.assertAlmostEqual(p,.05,places=12)
        a,_=benchmark_power(.2,2,7,39,.05)
        b,_=benchmark_power(.8,.5,7,39,.05)
        self.assertAlmostEqual(a,b,places=12)
        self.assertGreater(a,p)
        self.assertLess(benchmark_power(.2,2,7,39,3e-7)[0],a)

    def test_zero_rejection_interval_has_nonzero_upper_bound(self):
        low,high=wilson(0,40)
        self.assertAlmostEqual(low,0,places=12)
        self.assertGreater(high,0)
        self.assertEqual(wilson(0,0),(None,None))


if __name__=='__main__':unittest.main()
