"""Tests for dr4_bootstrap.py (standard library only). Apache-2.0."""
import io
import os
import tempfile
import unittest
import zipfile

import dr4_bootstrap as b


def make_zip(members):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, data in members.items():
            z.writestr(name, data)
    return buf.getvalue()


class TestBootstrap(unittest.TestCase):
    def test_archive_url(self):
        self.assertEqual(b.archive_url("v0.1.0"),
                         "https://github.com/munimtechnologies/proxima-c-dr4/archive/refs/tags/v0.1.0.zip")

    def test_extract_single_top(self):
        blob = make_zip({"proxima-c-dr4-0.1.0/proxima_job.py": "print(1)\n",
                         "proxima-c-dr4-0.1.0/gaiaorb/__init__.py": ""})
        with tempfile.TemporaryDirectory() as d:
            root = b.extract(blob, d)
            self.assertEqual(os.path.basename(root), "proxima-c-dr4-0.1.0")
            self.assertTrue(os.path.exists(os.path.join(root, "proxima_job.py")))
            self.assertTrue(os.path.exists(os.path.join(root, "gaiaorb", "__init__.py")))

    def test_rejects_traversal(self):
        with self.assertRaises(ValueError):
            b.safe_names(["top/../../evil.py"])
        with self.assertRaises(ValueError):
            b.safe_names(["/etc/passwd"])

    def test_rejects_two_tops(self):
        blob = make_zip({"a/x.py": "", "b/y.py": ""})
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                b.extract(blob, d)

    def test_split_args(self):
        args, rest = b.split_args(["--tag", "v0.1.0", "--", "--data-dir", "/data", "--validated"])
        self.assertEqual(args.tag, "v0.1.0")
        self.assertFalse(args.skip_tests)
        self.assertEqual(rest, ["--data-dir", "/data", "--validated"])
        args, rest = b.split_args([])
        self.assertEqual(args.tag, "v0.2.1")
        self.assertEqual(rest, [])


if __name__ == "__main__":
    unittest.main()
