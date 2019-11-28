import os, sys
sys.path.append("")

import Adiabatic_wind_solver as aws
import pytest
import numpy as np

def test_pc():
    '''Simple test of the percent change function'''
    s=aws.solver()
    assert s.percentChange(np.array([10,0]),np.array([11,0]))==10.0
    
def test_dw():
    '''Checks that the isothermal (base) case with gamma=1 has no change in temperature'''
    s=aws.solver()
    assert s.dw(1,[2,3,4])==0
    
def test_generateFunc():
    '''Tests that the RK integrator is working correctly for the isothermal case'''
    s=aws.solver()
    assert s.generateFunc(np.array([1,.001,0]))[2][-1]==pytest.approx(0.7145761164770619)

def test_findV0():
    '''Checks that the correct v0 is found for the adiabatic and isothermal cases'''
    s=aws.solver()
    assert s.findV0(.004,.006,.0001)==pytest.approx(0.005075319385528564)
    s1=aws.solver(5/3,.1)
    assert s1.findV0(.001,.003,.0001)==pytest.approx(0.0024991476535797122)
    
