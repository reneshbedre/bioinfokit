from bioinfokit.analys import norm, get_data, genfam
import pytest
import numpy as np
import pandas as pd
from unittest import TestCase
import urllib.request


class TestNormalization(TestCase):
    def test_cpm(self):
        df = get_data('sc_exp').data
        df = df.drop(['length'], axis=1)
        df = df.set_index('gene')
        nm = norm()
        nm.cpm(df=df)
        np.testing.assert_array_equal(round(nm.cpm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([100.695004, 101.731189, 74.721094, 92.633828, 74.270713, 95.314714])
                                      .round(2))

    def test_rpkm(self):
        df = get_data('sc_exp').data
        df = df.set_index('gene')
        nm = norm()
        nm.rpkm(df=df, gl='length')
        np.testing.assert_array_equal(round(nm.rpkm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([50.804745, 51.327542, 37.699846, 46.737552, 37.472610, 48.090169])
                                      .round(2))

    def test_tpm(self):
        df = get_data('sc_exp').data
        df = df.set_index('gene')
        nm = norm()
        nm.tpm(df=df, gl='length')
        np.testing.assert_array_equal(round(nm.tpm_norm.iloc[0], 2).to_numpy(),
                                      np.asarray([99.730156, 97.641941, 72.361658, 89.606265, 69.447237, 90.643338])
                                      .round(2))

    def test_genfam(self):
        id_data = pd.read_csv('https://reneshbedre.github.io/assets/posts/genfam/grai_id.txt', header=None)
        res = genfam()
        res.fam_enrich(id_file=id_data, species='grai', id_type=1)
        np.testing.assert_array_equal(res.df_enrich.iloc[0][1], 'MYB')
