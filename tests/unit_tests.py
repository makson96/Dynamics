from unittest import TestCase
import os
from pymol_plugin_dynamics import get_dynamics_dir, get_project_dirs, standardize_new_line_char


class Test(TestCase):
    def test_get_dynamics_dir(self):
        exp = "/home/{}/.dynamics/".format(os.getenv("USER"))
        obs = get_dynamics_dir()
        self.assertEqual(exp, obs)

    def test_get_project_dirs(self):
        value1 = "abcdefg"
        exp = "/home/{}/.dynamics/{}/".format(os.getenv("USER"), value1)
        obs = get_project_dirs(value1)
        self.assertEqual(exp, obs)

    def test_standardize_new_line_char1(self):
        exp = "aaa"
        obs = standardize_new_line_char("aaa")
        self.assertEqual(exp, obs)

    def test_standardize_new_line_char2(self):
        exp = "aaa\nbbb"
        obs = standardize_new_line_char("aaa\nbbb")
        self.assertEqual(exp, obs)

    def test_standardize_new_line_char3(self):
        exp = "aaa\nbbb\nccc"
        obs = standardize_new_line_char("aaa\rbbb\r\nccc")
        self.assertEqual(exp, obs)


if __name__ == "__main__":
    u_tests = Test()
    u_tests.test_get_dynamics_dir()
    u_tests.test_get_project_dirs()
    u_tests.test_standardize_new_line_char1()
    u_tests.test_standardize_new_line_char2()
    u_tests.test_standardize_new_line_char3()
