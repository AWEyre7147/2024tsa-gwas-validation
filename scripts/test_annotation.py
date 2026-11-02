"""Focused interval, orientation and checksum tests for the annotation workflow."""
from pathlib import Path
import hashlib
import tempfile
import unittest
from annotate_candidates import attributes, distance, flank_check, transcript_classification
from download_annotation_assets import published_md5


class Sequence:
    def __init__(self, value):
        self.value = value

    def get(self, chromosome, start, end):
        return self.value[start:end]


class AnnotationChecks(unittest.TestCase):
    def test_codon_changes_and_spliced_reverse_sequence(self):
        from annotate_coding_consequences import consequence, coding_sequence
        self.assertEqual(consequence('GAA', 'GAG'), ('E', 'E', 'synonymous_variant'))
        self.assertEqual(consequence('GAA', 'GCA'), ('E', 'A', 'missense_variant'))
        self.assertEqual(consequence('TGG', 'TGA'), ('W', '*', 'stop_gained'))
        self.assertEqual(coding_sequence(Sequence('ATGACCGT'), 'x', '+', [(1, 4, 0), (7, 8, 2)], 7), ('ATGAGT', 4))
        self.assertEqual(coding_sequence(Sequence('ACTCCCAT'), 'x', '-', [(1, 2, 2), (5, 8, 0)], 2), ('ATGGGT', 4))

    def test_closed_gene_boundaries(self):
        gene = {'start': 100, 'end': 200}
        self.assertEqual([distance(p, gene) for p in [99, 100, 200, 201]], [1, 0, 0, 1])

    def test_transcript_cds_utr_and_splice_boundaries(self):
        t = {'exon': [(100, 150), (200, 250)], 'CDS': [(120, 150), (200, 230)], 'UTR': []}
        expected = {100: 'UTR_overlap', 120: 'CDS_overlap', 150: 'CDS_overlap',
                    151: 'splice_boundary_2bp', 152: 'splice_boundary_2bp', 153: 'intron_overlap',
                    197: 'intron_overlap', 198: 'splice_boundary_2bp', 199: 'splice_boundary_2bp',
                    200: 'CDS_overlap', 250: 'UTR_overlap'}
        for pos, label in expected.items():
            self.assertEqual(transcript_classification(pos, t), label)
        self.assertEqual(transcript_classification(100, {'exon': [(100, 150)], 'CDS': [], 'UTR': []}),
                         'noncoding_transcript_exon')

    def test_reverse_complement_context_and_excluded_snp(self):
        from reference_alignment import COMP
        source = 'ACGTA' * 10 + 'C'
        target = source.translate(COMP)[::-1]
        target = target[:25] + 'A' + target[26:]
        config = {'context_flank_bp': 25, 'minimum_context_bases': 45, 'maximum_context_mismatches': 2}
        self.assertEqual(flank_check(Sequence(source), Sequence(target), 'x', 26, 'y', 26, '-', config), (50, 0, True))
        bad = 'N' * 51
        self.assertEqual(flank_check(Sequence(source), Sequence(bad), 'x', 26, 'y', 26, '-', config), (0, 0, False))

    def test_provider_attribute_escaping(self):
        self.assertEqual(attributes('ID=gene-x;Name=RNA%3Bname', 'RefSeq')['Name'], 'RNA;name')
        self.assertEqual(attributes('gene_id "ENSG1"; gene_biotype "lncRNA";', 'Ensembl')['gene_biotype'], 'lncRNA')

    def test_published_checksum_absolute_and_relative_paths(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'genome'
            path.write_bytes(b'genome')
            listing = Path(folder) / 'checksums'
            digest = hashlib.md5(b'genome').hexdigest()
            for name in ['./genome', '/mirrordata/hubs/genome']:
                listing.write_text(digest + '  ' + name + '\n')
                published_md5(path, listing, 'genome')
            path.write_bytes(b'changed')
            with self.assertRaises(ValueError):
                published_md5(path, listing, 'genome')


if __name__ == '__main__':
    unittest.main()
