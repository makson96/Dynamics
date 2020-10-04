from unittest import TestCase
from pymol_plugin_dynamics import get_dynamics_dir
import os


class Test(TestCase):
    def test_get_dynamics_dir(self):
        exp = "/home/{}/.dynamics/".format(os.getenv("USER"))
        obs = get_dynamics_dir()
        self.assertEqual(exp, obs)
