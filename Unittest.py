import unittest
from Particle import particle

class TestParticleMethods(unittest.TestCase):
    def test_negative_mass(self):
        with self.assertRaises(ValueError):
            Proton = particle(mass=-20)

    def test_invalid_position(self):
        with self.assertRaises(ValueError):
            Proton = particle(position='invalid_position')

    def test_invalid_velocity(self):
        with self.assertRaises(ValueError):
            Proton = particle(velocity='invalid_velocity')

    def test_invalid_force(self):
        with self.assertRaises(ValueError):
            Proton = particle(Force='invalid_force')

    def test_invalid_charge(self):
        with self.assertRaises(ValueError):
            Proton = particle(q=None)

    def test_invalid_name(self):
        with self.assertRaises(ValueError):
            Proton = particle(name=123)

    def test_invalid_B(self):
        with self.assertRaises(ValueError):
            Proton = particle(B='invalid_B')

    def test_negative_radius(self):
        with self.assertRaises(ValueError):
            Proton = particle(r_Cycl=-10)

    def test_negative_gap(self):
        with self.assertRaises(ValueError):
            Proton = particle(d=-10)

    def test_gap_greater_than_radius(self):
        with self.assertRaises(ValueError):
            Proton = particle(d=20, r_Cycl=10)

if __name__ == '__main__':
    unittest.main()
