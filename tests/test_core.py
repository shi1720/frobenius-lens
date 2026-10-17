import unittest,sys,random
from pathlib import Path
sys.path[:0]=[str(Path(__file__).resolve().parents[1]/'src'),str(Path(__file__).resolve().parent)]
from frobenius_lens import Experiment, certificate, recover, snapshot, minimum_repair
from independent import evolve

class TestLens(unittest.TestCase):
    def test_reference_example(self):
        es=[Experiment(n,64) for n in [32,64,128]]
        c=certificate(es,4)
        self.assertEqual((c['degree'],c['rank'],c['fiber_size']),(2,2,128))
        self.assertEqual(c['L'],[1,0,1])
        self.assertTrue(certificate(es+[Experiment(9,64)],4)['identifiable'])
        self.assertFalse(certificate(es+[Experiment(256,64)],4)['identifiable'])
    def test_empty_and_zero(self):
        self.assertEqual(certificate([],2)['fiber_size'],32)
        self.assertEqual(certificate([Experiment(9,64,(0,))],2)['rank'],0)
        self.assertEqual(recover([],[],1)['particular_rule'],[0,0,0])
    def test_constant(self):
        for n in [7,8,9,16]:
            self.assertEqual(certificate([Experiment(n,64,(1,)*n)],4)['rank'],1)
    def test_duplicate(self):
        e=Experiment(9,64)
        self.assertEqual(certificate([e],4)['L'],certificate([e,e],4)['L'])
    def test_recovery(self):
        rng=random.Random(47)
        for p in [2,3,5]:
            for r in range(5):
                for n in [2,5,7,11]:
                    for a in [0,1,2]:
                        es=[Experiment(n,p**a),Experiment(n+1,p**a)]
                        f=[rng.randrange(p) for _ in range(2*r+1)]
                        ys=[evolve(f,e.n,e.t,e.seed,r,p) for e in es]
                        z=recover(es,ys,r,p)
                        self.assertEqual([snapshot(z['particular_rule'],e,r,p) for e in es],ys)
                        if z['certificate']['identifiable']:
                            self.assertEqual(z['particular_rule'],f)
    def test_inconsistent(self):
        with self.assertRaises(ValueError):
            recover([Experiment(5,1,(0,))],[(1,0,0,0,0)],1)
    def test_observation_shape(self):
        with self.assertRaises(ValueError):
            recover([Experiment(5,1)],[(1,0)],1)
    def test_unsupported(self):
        for e in [Experiment(5,0),Experiment(5,3),Experiment(-1,1),Experiment(2,1,(1,1,1))]:
            with self.assertRaises(ValueError):
                certificate([e],1)
    def test_boundary(self):
        for r in range(1,5):
            self.assertFalse(certificate([Experiment(2*r,1)],r)['identifiable'])
            self.assertTrue(certificate([Experiment(2*r+1,1)],r)['identifiable'])
    def test_repair_no_solution(self):
        self.assertIsNone(minimum_repair([],4,2,[Experiment(8,64)]))

if __name__=='__main__':unittest.main()
