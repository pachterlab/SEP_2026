"""Regression checks for R source and argument handling across R launchers."""

import shutil
import subprocess
import unittest

from lib.common import run_r


@unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
class RRuntimeTests(unittest.TestCase):
    def test_source_escapes_and_spaced_arguments(self):
        code = r'''
        a <- commandArgs(TRUE)
        stopifnot(a[1] == "argument with spaces")
        cat(sub("factor\\(md\\$donor\\)", "", "factor(md$donor)7"), "\t", "ok", sep="")
        '''
        result = run_r(code, "argument with spaces", text=True, capture_output=True)
        self.assertEqual(result.stdout, "7\tok")

    def test_r_errors_propagate(self):
        with self.assertRaises(subprocess.CalledProcessError):
            run_r('stop("intentional test error")', text=True, capture_output=True)


if __name__ == "__main__":
    unittest.main()
