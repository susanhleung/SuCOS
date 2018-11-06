"""Tests for calc_SuCOS.py"""
import unittest
import os
import calc_SuCOS
from rdkit import Chem

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestCase(unittest.TestCase):
    
    def test1_SuCOS(self):
        """Check for perfect overlap with itself"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_lig.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score = float(ms[0].GetProp("SuCOS_score"))
        assert isclose(SuCOS_score, 1.0, abs_tol=1e-2)

    def test2_SuCOS(self):
        """Check Hydrogens don't affect score"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_ligH.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]
        
        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score = float(ms[0].GetProp("SuCOS_score"))
        assert isclose(SuCOS_score, 1.0, abs_tol=1e-2)

    def test3_SuCOS(self):
        """Test multiple prb molecules"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/mols.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 6
        for m in ms:
            SuCOS_score = float(m.GetProp("SuCOS_score"))
            assert SuCOS_score != 0

    def test5_SuCOS(self):
        """Prb molecule is same as ref but with extra group"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_lig_Ph.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score = float(ms[0].GetProp("SuCOS_score"))
        assert isclose(SuCOS_score, 1.0, abs_tol=1e-2)

    def test4_SuCOS(self):
        """Test rotation of COOH group. NB it lowers score by approx 0.05"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_lig_rotatedCOOH.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score = float(ms[0].GetProp("SuCOS_score"))
        assert isclose(SuCOS_score, 0.95, abs_tol=1e-2)

    def test6_SuCOS(self):
        """Complains if wrong input file format"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_lig.pdb"

        with self.assertRaises(ValueError): calc_SuCOS.main(ref_sdf, prb_sdf)


    def test7_SuCOS(self):
        """Prb molecule same as ref but with OH missing"""
        ref_sdf = "test_data/4e3g_lig.sdf"
        prb_sdf = "test_data/4e3g_lig_OHmissing.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score = float(ms[0].GetProp("SuCOS_score"))
        assert SuCOS_score < 1

        # More groups missing
        prb_sdf = "test_data/benzene.sdf"

        calc_SuCOS.main(ref_sdf, prb_sdf)
        outfile = "%s_SuCOS_score.sdf" % os.path.splitext(prb_sdf)[0]

        assert os.path.isfile(outfile)
        ms = Chem.SDMolSupplier(outfile)
        assert len(ms) == 1
        SuCOS_score2 = float(ms[0].GetProp("SuCOS_score"))
        assert SuCOS_score2 < 1
        assert SuCOS_score2 < SuCOS_score 

if __name__ == '__main__':
    print("Testing SuCOS")
    unittest.main()

