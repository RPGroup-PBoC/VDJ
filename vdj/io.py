"""
io.py

This submodule contains an array of utilities for simple file input/output and
other file parsing utilities 
"""
import numpy as np
import pandas as pd
import os

def endogenous_seqs():
    """
    Returns a dictionary of the sequence identity for . 
    """
    seqs = {'consensus': 'CACAGTGCTACAGACTGGAACAAAAACC',
            'WT1223rss': None,
            'DFL161':   None,
            'DFL1613': None,
            'V1-135': None,
            'V9-120': None,
            'V10-96': None,
            'V10-95': None,
            'V19-93': None}
    return seqs

def nucleotide_idx():
    """
    Returns the dictionary linking base identity to an integer
    """
    return {'A':0, 'C':1, 'G':2, 'T':3}

class ProcessTPM(object):
    R"""
    Base class for processing information from a compiled TPM experiment.
    """
    def __init__(self, file):
        R"""
        Parameters
        ---------
        file: str
            Path to the experimental `.mat` file. 
        """
        return None
