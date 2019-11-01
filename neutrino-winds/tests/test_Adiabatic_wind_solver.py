import os, sys
sys.path.append("")

import Adiabatic_wind_solver as aws
import pytest
import numpy as np

def test_pc():
    s=aws.solver()
    assert s.percentChange(np.array([10,0]),np.array([11,0]))==10.0
    
def test_dw():
    s=aws.solver()
    assert s.dw(1,[2,3,4])==0
    
def test_generateFunc():
    s=aws.solver()
    assert s.generateFunc(np.array([1,.001,0]))[2][-1]==pytest.approx(0.7145761164770619)

def test_findV0():
    s=aws.solver()
    assert s.findV0(.004,.006,.0001)==pytest.approx(0.005075319385528564)
    