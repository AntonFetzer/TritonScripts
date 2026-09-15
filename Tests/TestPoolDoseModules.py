import unittest
from unittest.mock import patch
from Dependencies.PoolDoseModules import poolDoseModules

class TestPoolDoseModules(unittest.TestCase):
    def item(self,dose,error,n,unit="MeV/g"):
        return {"doseLayer-1":dict(dose=dose,error=error,entries=n,unit=unit,**{"non-zeros":n})}
    @patch("Dependencies.PoolDoseModules.readDoseModules")
    def test_entry_weighted_units(self,reader):
        reader.side_effect=[self.item(2.,.2,100),self.item(4.,.4,300)]
        r=poolDoseModules(["a","b"])["doseLayer-1"]
        self.assertAlmostEqual(r["dose"],3.5)
        self.assertAlmostEqual(r["error"],(20**2+120**2)**.5/400)
        self.assertEqual(r["unit"],"MeV/g")
        self.assertEqual(r["entries"],400)
        self.assertFalse(r["inconsistent"])
        self.assertAlmostEqual(r["birge_ratio"], (150**2+150**2)**.5/(20**2+120**2)**.5)
    @patch("Dependencies.PoolDoseModules.readDoseModules")
    def test_zero_dose(self,reader):
        reader.side_effect=[self.item(0.,0.,100),self.item(0.,0.,100)]
        r=poolDoseModules(["a","b"])["doseLayer-1"]
        self.assertEqual(r["dose"],0.)
        self.assertEqual(r["birge_ratio"],0.)
    @patch("Dependencies.PoolDoseModules.readDoseModules")
    def test_unit_mismatch(self,reader):
        reader.side_effect=[self.item(2.,.2,100),self.item(2.,.2,100,"rad")]
        with self.assertRaises(ValueError):poolDoseModules(["a","b"])
