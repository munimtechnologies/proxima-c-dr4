"""Tests for io_gaia.py v0.2.0 (DataLink data-structure rename; rejecting error responses).
Standard library only; no network. MV2 / Munim, Inc. Apache-2.0. Made with AI agents."""
import io
import os
import tempfile
import unittest
import zipfile

import io_gaia

HEADER = "source_id,transit_id,obs_time_tcb,obs_time_bary_corr,centroid_pos_al,centroid_pos_error_al," \
         "scan_pos_angle,parallax_factor_al,used_by_agis_al\n"
GOOD_CSV = (HEADER +
            '123,1,"(2.5e17, 2.5e17, 2.5e17)",1.0e8,"(0.1, 0.2, 0.3)","(0.5, 0.5, 0.5)",'
            '"(10.0, 10.0, 10.0)",0.3,"(true, true, true)"\n'
            '123,2,"(2.6e17, 2.6e17, 2.6e17)",-1.0e8,"(0.4, 0.5, 0.6)","(0.5, 0.5, 0.5)",'
            '"(70.0, 70.0, 70.0)",-0.2,"(true, true, true)"\n').encode()


def zip_blob(members):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in members.items():
            z.writestr(name, data)
    return buf.getvalue()


def quiet(_msg):
    pass


class TestUrls(unittest.TestCase):
    def test_new_structure_first_then_old_then_none(self):
        urls = io_gaia.datalink_urls(123)
        self.assertEqual(len(urls), 12)
        for u in urls:
            self.assertTrue(u.startswith("https://gea.esac.esa.int/data-server/data?"))
            self.assertIn("RETRIEVAL_TYPE=EPOCH_ASTROMETRY", u)
            self.assertIn("RELEASE=Gaia+DR4", u)
        for u in urls[:4]:
            self.assertIn("DATA_STRUCTURE=DATAMODEL_STANDARD", u)
        for u in urls[4:8]:
            self.assertIn("DATA_STRUCTURE=INDIVIDUAL", u)
        for u in urls[8:]:
            self.assertNotIn("DATA_STRUCTURE", u)
        self.assertIn("ID=Gaia+DR4+123", urls[0])
        self.assertIn("FORMAT=CSV", urls[0])
        self.assertIn("FORMAT=votable", urls[1])
        self.assertIn("ID=123&", urls[2])


class TestSaveEpochs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_good_csv_kept_and_readable(self):
        paths = io_gaia.save_epochs(GOOD_CSV, self.d, "x", log=quiet)
        self.assertEqual(len(paths), 1)
        rows = io_gaia.read_table_rows(paths[0])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["centroid_pos_al"], [0.1, 0.2, 0.3])
        self.assertEqual(rows[1]["used_by_agis_al"], [True, True, True])

    def test_error_text_rejected(self):
        blob = b"ERROR: Data structure INDIVIDUAL is not supported. Use DATAMODEL_STANDARD.\n" * 4
        self.assertEqual(io_gaia.save_epochs(blob, self.d, "x", log=quiet), [])

    def test_header_only_rejected(self):
        self.assertEqual(io_gaia.save_epochs(HEADER.encode(), self.d, "x", log=quiet), [])

    def test_html_rejected(self):
        blob = b"<!DOCTYPE html><html><body>Service Unavailable</body></html>" + b" " * 300
        self.assertEqual(io_gaia.save_epochs(blob, self.d, "x", log=quiet), [])

    def test_votable_error_rejected(self):
        blob = (b'<?xml version="1.0" encoding="UTF-8"?>\n<VOTABLE version="1.4"><RESOURCE type="results">'
                b'<INFO name="QUERY_STATUS" value="ERROR">Unknown data structure</INFO></RESOURCE></VOTABLE>')
        self.assertEqual(io_gaia.save_epochs(blob, self.d, "x", log=quiet), [])

    def test_zip_keeps_only_epoch_tables(self):
        blob = zip_blob({"EPOCH_ASTROMETRY/EPOCH_ASTROMETRY-Gaia DR4 123.csv": GOOD_CSV,
                         "README.txt": b"Gaia DataLink bundle",
                         "other.csv": b"a,b\n1,2\n"})
        paths = io_gaia.save_epochs(blob, self.d, "x", log=quiet)
        self.assertEqual([os.path.basename(p) for p in paths], ["EPOCH_ASTROMETRY-Gaia DR4 123.csv"])

    def test_zip_without_epochs_rejected(self):
        blob = zip_blob({"x.csv": b"a,b\n1,2\n"})
        self.assertEqual(io_gaia.save_epochs(blob, self.d, "x", log=quiet), [])

    def test_bad_zip_rejected(self):
        self.assertEqual(io_gaia.save_epochs(b"PK" + b"\x00" * 400, self.d, "x", log=quiet), [])


if __name__ == "__main__":
    unittest.main()
